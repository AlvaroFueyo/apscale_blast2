"""Deterministic assignment of one query's hits, independent of I/O."""

from collections import Counter

from .filtering import choose_flag_rest, thresholds_to_dict, trim_by_similarity

TAX_COLS = ["Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species"]
PLURALS = ["kingdoms", "phyla", "classes", "orders", "families", "genera", "species"]


def metric_key(row):
    return (-row["Similarity"], row["evalue"], -row["query_coverage"], tuple(row.get(c, "") for c in TAX_COLS))


def candidate_labels(rows):
    return sorted({next((r[c] for c in reversed(TAX_COLS) if r[c]), "") for r in rows} - {""})


def assign_query(qid, hits, opts):
    """Preserve rank compatibility without confusing an unknown rank with agreement."""
    base = {"unique ID": qid, **dict.fromkeys(TAX_COLS, "NoMatch"), "Similarity": 0.0,
            "query_coverage": 0.0, "evalue": 1.0, "Flag": "", "Ambiguous taxa": "",
            "assignment_status": "no_match", "assigned_rank": "", "hit_count": len(hits),
            "candidate_taxa": 0, "hit_limit_reached": len({r["Subject ID"] for r in hits}) >= opts.max_target_seqs}
    base["reference_uncertainty"] = ""
    base["reference_taxonomy_conflict"] = ""
    base["reference_missing_nodes"] = ""
    if not hits:
        return base
    selected = hits
    if opts.prefer_qcov > 0:
        preferred = [r for r in hits if r["query_coverage"] >= opts.prefer_qcov]
        selected = preferred or hits
    if opts.filter_mode == 1:
        similarity = max(r["Similarity"] for r in selected)
        selected = [r for r in selected if abs(r["Similarity"] - similarity) <= 1e-9]
        evalue = min(r["evalue"] for r in selected)
        selected = [r for r in selected if r["evalue"] == evalue]
    else:
        evalue = min(r["evalue"] for r in selected)
        selected = [r for r in selected if r["evalue"] == evalue]
        similarity = max(r["Similarity"] for r in selected)
        selected = [r for r in selected if abs(r["Similarity"] - similarity) <= 1e-9]
    thresholds = thresholds_to_dict(opts.thresholds)
    base["reference_uncertainty"] = ";".join(sorted({r.get("species_uncertainty", "") for r in selected} - {""}))
    base["reference_taxonomy_conflict"] = ";".join(sorted({r.get("reference_taxonomy_conflict", "") for r in selected} - {""}))
    base["reference_missing_nodes"] = ";".join(sorted({r.get("reference_missing_nodes", "") for r in selected} - {""}))
    trimmed = []
    for row in selected:
        row = row.copy()
        trim_by_similarity(row, similarity, thresholds)
        trimmed.append(row)
    unique = {}
    for row in sorted(trimmed, key=metric_key):
        unique.setdefault(tuple(row[c] for c in TAX_COLS), row)
    taxa = list(unique.values())
    out = {**base, **{c: taxa[0][c] for c in TAX_COLS + ["Similarity", "query_coverage", "evalue"]}}
    out["candidate_taxa"] = len(taxa)
    missing = any(not r.get("taxonomy_mapped", True) for r in selected)
    if opts.flag_scheme == "apscale" and len(taxa) > 1:
        if similarity < thresholds["Species"]:
            for rank in reversed(TAX_COLS):
                for row in taxa:
                    row[rank] = ""
                if len({tuple(r[c] for c in TAX_COLS) for r in taxa}) == 1:
                    break
            out.update({c: taxa[0][c] for c in TAX_COLS})
        else:
            counts = Counter(tuple(r[c] for c in TAX_COLS) for r in trimmed)
            dominant = [key for key, count in counts.items() if count == max(counts.values())]
            if len(dominant) == 1:
                out.update(dict(zip(TAX_COLS, dominant[0])))
                out["Flag"] = "F1 (Dominant species)"
                out["Ambiguous taxa"] = ", ".join(candidate_labels(taxa))
            else:
                chosen = choose_flag_rest(taxa, [r["Species"] for r in taxa if r["Species"]])
                out.update({c: chosen[c] for c in TAX_COLS + ["Flag", "Ambiguous taxa"]})
    elif opts.flag_scheme == "apscale2":
        consensus = dict.fromkeys(TAX_COLS, "")
        conflict = None
        for index, rank in enumerate(TAX_COLS):
            values = {r[rank] for r in taxa if r[rank]}
            if len(values) > 1:
                conflict = index
                break
            if values:
                consensus[rank] = next(iter(values))
        out.update(consensus)
        if conflict is not None:
            out["Ambiguous taxa"] = ", ".join(candidate_labels(taxa))
            if TAX_COLS[conflict] == "Species" and consensus["Genus"]:
                genus = consensus["Genus"]
                species = sorted({r["Species"] for r in taxa if r["Species"]})
                epithets = [s[len(genus) + 1:] if s.startswith(genus + " ") else s for s in species]
                out["Species"] = f"{genus} {'/'.join(epithets)}" if len(species) == 2 else f"{genus} sp."
                out["Flag"] = "Fl1 Two species of one genus" if len(species) == 2 else "Fl1 More than two species of one genus"
            else:
                out["Flag"] = f"Fl{7 - conflict} Two or more {PLURALS[conflict]} (trimming to MRCA)"
    if missing:
        out["assignment_status"] = "taxonomy_missing"
    elif out["Flag"] or out["reference_taxonomy_conflict"]:
        out["assignment_status"] = "ambiguous"
    elif not any(out[c] for c in TAX_COLS):
        out["assignment_status"] = "unresolved"
    else:
        out["assignment_status"] = "assigned"
    ranks = [c for c in TAX_COLS if out[c] and out[c] != "NoMatch"]
    if out["Flag"].startswith(("Fl1", "F2", "F3")) and "Species" in ranks:
        ranks.remove("Species")
    out["assigned_rank"] = ranks[-1] if ranks else ""
    return out
