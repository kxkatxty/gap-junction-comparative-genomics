"""Shared paths for the pannexin workspace (independent of the gap-junction thesis)."""

from __future__ import annotations

import csv
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = PROJECT_ROOT / "project" / "metadata"
DATA_DIR = PROJECT_ROOT / "project" / "data"
RESULTS_DIR = PROJECT_ROOT / "project" / "results"
REFERENCES_DIR = DATA_DIR / "references" / "pannexins"
ANNOTATIONS_DIR = DATA_DIR / "annotations" / "pannexin"
LOGS_DIR = PROJECT_ROOT / "project" / "logs"

# Chordate / vertebrate paralogs used as working labels.
PANNEXIN_TYPES = ("PANX1", "PANX2", "PANX3")

# Coarse clade palette for later figures (not purple-on-white).
CLADE_COLORS = {
    "Mammal": "#B45309",
    "Bird": "#0369A1",
    "Reptile": "#15803D",
    "Amphibian": "#0F766E",
    "Ray-finned fish": "#1D4ED8",
    "Lobe-finned fish": "#7C3AED",
    "Cartilaginous fish": "#0E7490",
    "Jawless vertebrate": "#64748B",
    "Tunicate / lancelet": "#78716C",
    "Other": "#94A3B8",
}

TYPE_COLORS = {
    "PANX1": "#C2410C",
    "PANX2": "#1D4ED8",
    "PANX3": "#15803D",
    "other/unknown": "#94A3B8",
}


def read_csv_rows(path: Path) -> list[dict[str, str]]:
    if not path.exists() or path.stat().st_size == 0:
        return []
    with path.open("r", encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def write_csv(path: Path, rows: list[dict[str, str]], fieldnames: list[str] | None = None) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fields = fieldnames or list(rows[0].keys())
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def classify_pannexin_type(text: str) -> str:
    t = (text or "").upper().replace("-", " ").replace("_", " ")
    t_compact = t.replace(" ", "")
    if "PANX1" in t_compact or "PX1" in t_compact or "PANNEXIN1" in t_compact:
        return "PANX1"
    if "PANX2" in t_compact or "PX2" in t_compact or "PANNEXIN2" in t_compact:
        return "PANX2"
    if "PANX3" in t_compact or "PX3" in t_compact or "PANNEXIN3" in t_compact:
        return "PANX3"
    # UniProt descriptions like "Pannexin 1a"
    if "PANNEXIN 1" in t or "PANNEXIN1" in t_compact:
        return "PANX1"
    if "PANNEXIN 2" in t or "PANNEXIN2" in t_compact:
        return "PANX2"
    if "PANNEXIN 3" in t or "PANNEXIN3" in t_compact:
        return "PANX3"
    return "other/unknown"
