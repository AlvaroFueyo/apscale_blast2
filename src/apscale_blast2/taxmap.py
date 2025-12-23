"""Taxonomy mapping loader.

Loads the taxonomy mapping table associated with a local BLAST database and
returns it as a dictionary for fast lookups.
"""

from __future__ import annotations
import os, glob, pandas as pd, re
from functools import lru_cache
from typing import List, Dict
try:
    import pyarrow.parquet as pq
except Exception:
    pq = None
from .taxonomy_clean import clean_genus, clean_species

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
    if not isinstance(val, str): return "" if pd.isna(val) else str(val)
    v = val.strip()
    v = _PREFIX_RE.sub("", v)
    return v

def _candidate_dirs(db_prefix: str):
    paths = []
    if os.path.isdir(db_prefix):
        paths.append(os.path.abspath(db_prefix))
        base_dir = os.path.abspath(db_prefix)
    else:
        base_dir = os.path.dirname(os.path.abspath(db_prefix))
        paths.append(base_dir)
    parent = os.path.dirname(base_dir); grand = os.path.dirname(parent)
    for d in [parent, grand]:
        if d and d not in paths: paths.append(d)
    return paths

def find_taxmaps_paths(db_prefix: str):
    paths = []
    for d in _candidate_dirs(db_prefix):
        for pat in ["db_taxonomy*.parquet*", "*taxonomy*.parquet*", "db_taxonomy*.csv*", "*taxonomy*.csv*"]:
            paths.extend(glob.glob(os.path.join(d, pat)))
    base = _candidate_dirs(db_prefix); pref = {d:i for i,d in enumerate(base)}
    paths = sorted(set(paths), key=lambda p: (0 if ".parquet" in os.path.basename(p).lower() else 1, pref.get(os.path.dirname(os.path.abspath(p)), 99), p))
    return paths

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
    keep = ["Sequence ID"] + CANON_RANKS
    out = out[keep].copy()
    for c in CANON_RANKS:
        out[c] = out[c].map(_clean_rank)
    out["genus_name"] = out["genus_name"].map(clean_genus)
    out["species_name"] = out.apply(lambda r: clean_species(r.get("species_name","")) or "", axis=1)
    return out

@lru_cache(maxsize=32)
def load_taxmap_as_dict(db_prefix: str):
    last_err = None
    for p in find_taxmaps_paths(db_prefix):
        try:
            if p.lower().endswith(".parquet") or ".parquet" in os.path.basename(p).lower():
                if pq is None: 
                    last_err = "pyarrow not available"; continue
                df = pq.read_table(p).to_pandas()
            else:
                df = pd.read_csv(p)
            df = _normalize_columns(df)
            return {str(r["Sequence ID"]): [str(r[c]) for c in ["kingdom_name","phylum_name","class_name","order_name","family_name","genus_name","species_name"]] for _, r in df.iterrows()}
        except Exception as e:
            last_err = str(e); continue
    raise FileNotFoundError(f"No taxonomy parquet/csv found near the database ({db_prefix}). Last error: {last_err}")
