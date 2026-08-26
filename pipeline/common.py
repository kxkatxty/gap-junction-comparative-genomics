"""Shared paths and helpers for the gap-junction pipeline."""

from __future__ import annotations

import csv
import json
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
METADATA_DIR = PROJECT_ROOT / "project" / "metadata"
DATA_DIR = PROJECT_ROOT / "project" / "data"
RESULTS_DIR = PROJECT_ROOT / "project" / "results"
REFERENCES_DIR = DATA_DIR / "references"
ANNOTATIONS_DIR = DATA_DIR / "annotations"
SYNVOY_DIR = PROJECT_ROOT / "SynVoy"

DISCOVERY_DIRS = (
    RESULTS_DIR / "innexin_discovery",
    RESULTS_DIR / "connexin_discovery",
    RESULTS_DIR / "connexin_evolutionary_discovery",
    RESULTS_DIR / "connexin_evolutionary_batch2_discovery",
)

FAMILIES = ("innexin", "connexin", "pannexin")
REJECTED_RANK = "rejected_false_positive"


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


def infer_family_from_path(path: Path) -> str:
    text = str(path).lower()
    if "innexin" in text:
        return "innexin"
    if "connexin" in text:
        return "connexin"
    return "unknown"


def iter_discovery_summaries() -> list[tuple[Path, dict, str]]:
    rows: list[tuple[Path, dict, str]] = []
    for root in DISCOVERY_DIRS:
        if not root.exists():
            continue
        family = infer_family_from_path(root)
        for summary_path in sorted(root.glob("*/discovery_summary.json")):
            data = json.loads(summary_path.read_text(encoding="utf-8"))
            rows.append((summary_path, data, family))
    return rows
