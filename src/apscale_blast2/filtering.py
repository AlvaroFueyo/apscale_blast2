"""Filtering and flagging helpers.

Helpers to apply similarity-threshold trimming and the APSCALE-like F1–F4 flag logic.
"""

from __future__ import annotations
from typing import Dict, List

def thresholds_to_dict(thr_str: str | None) -> Dict[str,int]:
    defaults = ['97','95','90','87','85']
    parts = [t.strip() for t in (thr_str or ','.join(defaults)).split(',')]
    if len(parts) != 5: parts = defaults
    s,g,f,o,c = [int(x) for x in parts]
    return {"Species": s, "Genus": g, "Family": f, "Order": o, "Class": c}

def trim_by_similarity(row, max_sim: float, thr: Dict[str,int]):
    if max_sim >= thr["Species"]:
        return
    if max_sim < thr["Species"] and max_sim >= thr["Genus"]:
        row["Species"] = ""; return
    if max_sim < thr["Genus"] and max_sim >= thr["Family"]:
        row["Species"] = ""; row["Genus"] = ""; return
    if max_sim < thr["Family"] and max_sim >= thr["Order"]:
        row["Species"] = ""; row["Genus"] = ""; row["Family"] = ""; return
    if max_sim < thr["Order"] and max_sim >= thr["Class"]:
        row["Species"] = ""; row["Genus"] = ""; row["Family"] = ""; row["Order"] = ""; return
    for k in ["Species","Genus","Family","Order","Class"]:
        row[k] = ""

def choose_flag_rest(rows_dedup, ambiguous_species):
    if len(rows_dedup)==1:
        r=rows_dedup[0]; r["Flag"]=""; r["Ambiguous taxa"] = ""
        return r
    genera=set([r["Genus"] for r in rows_dedup if r["Genus"]])
    species=sorted(set([r for r in ambiguous_species if r]))
    if len(genera)==1:
        if len(species)==2:
            genus=list(genera)[0]; spec="/".join([s.replace(genus+" ","") for s in species])
            r=rows_dedup[0].copy(); r["Species"]=f"{genus} {spec}"; r["Flag"]="F2 (Two species of one genus)"; r["Ambiguous taxa"]=", ".join(species); return r
        if len(species)>2:
            r=rows_dedup[0].copy(); r["Species"]=f"{list(genera)[0]} sp."; r["Flag"]="F3 (Multiple species of one genus)"; r["Ambiguous taxa"]=", ".join(species); return r
    levels = ["Species","Genus","Family","Order","Class","Phylum","Kingdom"]
    tmp=[r.copy() for r in rows_dedup]
    for lv in levels:
        for r in tmp: r[lv]=""
        seen=set(); ded=[]
        for r in tmp:
            key=(r["Kingdom"],r["Phylum"],r["Class"],r["Order"],r["Family"],r["Genus"],r["Species"])
            if key not in seen: seen.add(key); ded.append(r)
        tmp=ded
        if len(tmp)==1:
            r=tmp[0]; r["Flag"]="F4 (Trimming to MRCA)"; r["Ambiguous taxa"]=", ".join(species); return r
    r=tmp[0]; r["Flag"]="F4 (Trimming to MRCA)"; r["Ambiguous taxa"]=", ".join(species); return r
