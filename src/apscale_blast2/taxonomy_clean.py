"""Taxonomy string normalization.

Helpers to normalise taxonomy names before writing outputs.
"""

from __future__ import annotations
import re

RE_CANDIDATUS = re.compile(r"^\s*(?:candidatus)\s+", re.I)
RE_QUAL = re.compile(r"(?:sp|spp|cf|aff|nr|complex|group|uncultured|unverified|environmental|metagenome|metagenomic|bacterium|archaeon|eukaryote)\.?,?", re.I)
RE_MULTI = re.compile(r"[,/;]+")
RE_GENUS = re.compile(r"^[A-Z][a-zA-Z-]+$")
RE_BINOMIAL = re.compile(r"^\s*([A-Z][a-zA-Z-]+)\s+([a-z][a-zA-Z-]+)")
RE_PREFIX = re.compile(r"^(?:kingdom|superkingdom|phylum|class|order|family|genus|species)\s+", re.I)

_PLACEHOLDER_VALUES = {
    "", "na", "n/a", "none", "null", "nan",
    "unclassified", "unidentified", "unknown", "uncultured",
    "unassigned", "unresolved", "environmental sample", "environmental samples",
    "metagenome", "metagenomic", "incertae sedis", "other", "miscellaneous",
    "bacterium", "archaeon", "eukaryote",
    "kingdom", "phylum", "class", "order", "family", "genus", "species",
}

def clean_taxon_name(v: str) -> str:
    """Generic cleaner for higher-rank taxonomy labels."""
    if not isinstance(v, str):
        return ""
    s = v.strip()
    if not s:
        return ""
    s = RE_PREFIX.sub("", s).strip()
    s = RE_CANDIDATUS.sub("", s).strip()
    if not s:
        return ""
    if s.lower() in _PLACEHOLDER_VALUES:
        return ""
    return s


def clean_genus(g: str) -> str:
    if not isinstance(g, str):
        return ""
    g = clean_taxon_name(g)
    if not g:
        return ""
    g = RE_QUAL.sub("", g).strip()
    return g if RE_GENUS.match(g) else ""


def clean_species(s: str) -> str:
    if not isinstance(s, str):
        return ""
    s = clean_taxon_name(s)
    if not s:
        return ""
    if RE_MULTI.search(s):
        return ""
    s = RE_QUAL.sub("", s).strip()
    m = RE_BINOMIAL.match(s)
    if not m:
        return ""
    genus, epithet = m.group(1), m.group(2)
    if not RE_GENUS.match(genus):
        return ""
    if RE_QUAL.fullmatch(epithet):
        return ""
    return f"{genus} {epithet}"
