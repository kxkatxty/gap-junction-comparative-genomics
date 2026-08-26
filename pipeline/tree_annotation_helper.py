#!/usr/bin/env python3
"""Prepare tree tip labels with family/species/clade metadata for plotting."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from pathlib import Path

from pipeline.common import METADATA_DIR, write_csv

DEFAULT_OUT = METADATA_DIR / "tree_tip_labels.csv"


def parse_newick_tips(tree_path: Path) -> list[str]:
    text = tree_path.read_text(encoding="utf-8").strip()
    # strip branch lengths; keep tip names before ':'
    tips = re.findall(r"([A-Za-z0-9_.|:-]+):[0-9.eE+-]+", text)
    if not tips:
        tips = [part.split(":")[0] for part in text.split(",") if part and not part.startswith("(")]
    return [t.strip("();") for t in tips if t.strip("();")]


def load_metadata() -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    master = METADATA_DIR / "thesis_master_summary.csv"
    if not master.exists():
        return lookup
    with master.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle):
            lookup[row["organism"]] = row
    return lookup


def build_labels(tree_path: Path) -> list[dict[str, str]]:
    meta = load_metadata()
    rows: list[dict[str, str]] = []
    for tip in parse_newick_tips(tree_path):
        organism = tip.split("|", 1)[0]
        info = meta.get(organism, {})
        rows.append(
            {
                "tip_label": tip,
                "organism": organism,
                "family": info.get("family", ""),
                "panel_role": info.get("panel_role", ""),
                "annotation_category": info.get("annotation_category", ""),
                "itol_label": f"{organism} [{info.get('family', '?')}]",
            }
        )
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Annotate phylogenetic tree tip labels.")
    parser.add_argument("tree", type=Path, help="Newick tree file.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.tree.exists():
        print(f"Tree not found: {args.tree}", file=sys.stderr)
        return 1
    rows = build_labels(args.tree)
    write_csv(args.output, rows)
    print(f"Wrote {len(rows)} tip labels -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
