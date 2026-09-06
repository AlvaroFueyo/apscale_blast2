"""Taxonomy mapping loader.

Loads the taxonomy mapping table associated with a local BLAST database and
returns it as a dictionary for fast lookups.
"""

from __future__ import annotations
import os, glob, pandas as pd, re
try:
    import pyarrow.parquet as pq
except Exception:
    pq = None
from .taxonomy_clean import clean_genus, clean_species, clean_taxon_name, species_uncertainty

CANON_RANKS = ["kingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]
RANK_SYNONYMS = {
    "kingdom_name": ["kingdom_name","Kingdom","kingdom","domain","superkingdom","superkingdom_name"],
    "phylum_name":  ["phylum_name","Phylum","phylum","division","division_name"],
    "class_name":   ["class_name","Class","class"],
    "order_name":   ["order_name","Order","order"],
    "family_name":  ["family_name","Family","family"],
    "genus_name":   ["genus_name","Genus","genus"],
    "species_name": ["species_name","Species","species","species_name_corrected"],
}
ID_CANDS = ["Sequence ID","sequence_id","seqid","id","header","Accession","accession"]
_PREFIX_RE = re.compile(r'^(kingdom|superkingdom|phylum|class|order|family|genus|species)\s+', re.I)

def _clean_rank(val: str) -> str:
    if not isinstance(val, str):
        return "" if pd.isna(val) else str(val)
    v = val.strip()
    v = _PREFIX_RE.sub("", v)
    return clean_taxon_name(v)

def _candidate_dirs(db_prefix: str):
    base = os.path.abspath(db_prefix)
    if not os.path.isdir(base):
        base = os.path.dirname(base)
    paths = [base]
    # Builders put indices in <database>/db and the taxonomy beside that folder.
    if os.path.basename(base).lower() == "db":
        paths.append(os.path.dirname(base))
    return paths

def find_taxmaps_paths(db_prefix: str):
    for d in _candidate_dirs(db_prefix):
        paths = set()
        for pat in ["*taxonomy*.parquet*", "*taxonomy*.csv", "*taxonomy*.csv.gz", "*taxonomy*.tsv", "*taxonomy*.tsv.gz"]:
            paths.update(glob.glob(os.path.join(d, pat)))
        if paths:
            canonical = [p for p in paths if os.path.basename(p).startswith("db_taxonomy.")]
            candidates = canonical or sorted(paths)
            if len(candidates) != 1:
                raise ValueError(f"Ambiguous taxonomy tables in {d}: {sorted(candidates)}")
            return candidates
    return []

def _normalize_columns(df: pd.DataFrame) -> pd.DataFrame:
    cols = list(df.columns); low = {c.lower(): c for c in cols}
    id_col = None
    for c in ID_CANDS:
        if c in df.columns: id_col = c; break
        if c.lower() in low: id_col = low[c.lower()]; break
    if id_col is None: raise ValueError("No 'Sequence ID' column found in taxonomy table.")
    out = df.rename(columns={id_col: "Sequence ID"}).copy()
    for dst, alts in RANK_SYNONYMS.items():
        found = None
        for a in alts:
            if a in out.columns: found = a; break
            if a.lower() in low: found = low[a.lower()]; break
        if found is None:
            out[dst] = ""
        else:
            out.rename(columns={found: dst}, inplace=True)
    original = out["species_original"].fillna("").astype(str) if "species_original" in out else out["species_name"].fillna("").astype(str)
    conflicts = out["taxonomy_conflict"].fillna("").astype(str).copy() if "taxonomy_conflict" in out else pd.Series("", index=out.index)
    missing_nodes = out["taxonomy_missing_nodes"].fillna("").astype(str) if "taxonomy_missing_nodes" in out else ""
    keep = ["Sequence ID"] + CANON_RANKS
    out = out[keep].copy()
    out["species_original"] = original
    out["species_uncertainty"] = out["species_original"].map(species_uncertainty)
    for c in CANON_RANKS:
        out[c] = out[c].map(_clean_rank)
    out["genus_name"] = out["genus_name"].map(clean_genus)
    out.loc[out["species_uncertainty"] == "hybrid", "genus_name"] = ""
    out["species_name"] = [clean_species(s, genus=g) for s, g in zip(out["species_name"], out["genus_name"])]
    mismatch = out["species_name"].ne("") & out["genus_name"].ne("") & out["species_name"].str.split().str[0].ne(out["genus_name"])
    out.loc[mismatch, "species_name"] = ""
    conflicts.loc[mismatch] = conflicts.loc[mismatch].map(lambda value: (value + "; " if value else "") + "species/genus mismatch")
    out["taxonomy_conflict"] = conflicts
    out["taxonomy_missing_nodes"] = missing_nodes
    return out

def load_taxmap_as_dict(db_prefix: str):
    paths = find_taxmaps_paths(db_prefix)
    if not paths:
        raise FileNotFoundError(f"No taxonomy table associated with database: {db_prefix}")
    path = paths[0]
    if ".parquet" in os.path.basename(path).lower():
        df = pd.read_parquet(path)
    else:
        df = pd.read_csv(path, sep="\t" if ".tsv" in path.lower() else ",", dtype=str, keep_default_na=False)
    df = _normalize_columns(df).fillna("")
    if df.empty or not df[CANON_RANKS].ne("").any(axis=None):
        raise ValueError(f"Taxonomy table has no informative ranks: {path}")
    result = {}
    for identifier, *ranks in df.itertuples(index=False, name=None):
        key = str(identifier).strip()
        if not key:
            raise ValueError(f"Empty reference ID in {path}")
        if key in result and result[key] != ranks:
            raise ValueError(f"Conflicting taxonomy for reference ID {key!r}: {path}")
        result[key] = ranks
    return result
