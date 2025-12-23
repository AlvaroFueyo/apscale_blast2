"""Taxonomy string normalization.

Small helpers to normalise genus/species names before writing outputs.
"""

from __future__ import annotations
import re
RE_CANDIDATUS = re.compile(r"^\s*(?:candidatus)\s+", re.I)
RE_QUAL = re.compile(r"\b(?:sp|spp|cf|aff|nr|complex|group|uncultured|unverified|environmental|metagenome|metagenomic|bacterium|archaeon|eukaryote)\b\.?", re.I)
RE_MULTI = re.compile(r"[,/;]+")
RE_GENUS = re.compile(r"^[A-Z][a-zA-Z-]+$")
RE_BINOMIAL = re.compile(r"^\s*([A-Z][a-zA-Z-]+)\s+([a-z][a-zA-Z-]+)\b")
def clean_genus(g: str) -> str:
    if not isinstance(g, str): return ""
    g = g.strip()
    g = RE_CANDIDATUS.sub("", g).strip()
    g = RE_QUAL.sub("", g).strip()
    return g if RE_GENUS.match(g) else ""
def clean_species(s: str) -> str:
    if not isinstance(s, str): return ""
    s = s.strip()
    s = RE_CANDIDATUS.sub("", s).strip()
    if RE_MULTI.search(s): return ""
    s = RE_QUAL.sub("", s).strip()
    m = RE_BINOMIAL.match(s)
    if not m: return ""
    genus, epithet = m.group(1), m.group(2)
    if not RE_GENUS.match(genus): return ""
    if RE_QUAL.fullmatch(epithet): return ""
    return f"{genus} {epithet}"
