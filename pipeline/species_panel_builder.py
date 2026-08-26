#!/usr/bin/env python3
"""Build a unified species panel table (reference vs discovery)."""

from __future__ import annotations

import argparse
import csv
import sys
from pathlib import Path

from pipeline.common import METADATA_DIR, PROJECT_ROOT, read_csv_rows, write_csv

DEFAULT_OUT = METADATA_DIR / "species_panel_master.csv"

PANEL_SOURCES = (
    ("reference_innexin", METADATA_DIR / "species_config.csv", "innexin"),
    ("reference_connexin", METADATA_DIR / "species_config.csv", "connexin"),
    ("discovery_innexin", METADATA_DIR / "not_annotated_innexin.txt", "innexin"),
    ("discovery_innexin_clade", METADATA_DIR / "species_discovery_queue.txt", "innexin"),
    ("discovery_connexin_evolutionary", METADATA_DIR / "not_annotated_connexin_evolutionary.txt", "connexin"),
    ("discovery_connexin_evolutionary_batch2", METADATA_DIR / "not_annotated_connexin_evolutionary_batch2.txt", "connexin"),
)


def load_species_txt(path: Path) -> list[str]:
    if not path.exists():
        return []
    return [
        line.strip()
        for line in path.read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.startswith("#")
    ]


def build_panel() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    seen: set[tuple[str, str, str]] = set()

    for panel_role, path, family in PANEL_SOURCES:
        if path.suffix == ".csv":
            for row in read_csv_rows(path):
                if row.get("family") != family and panel_role.startswith("reference"):
                    continue
                organism = row.get("organism", "").strip()
                if not organism:
                    continue
                key = (family, organism, panel_role)
                if key in seen:
                    continue
                seen.add(key)
                rows.append(
                    {
                        "family": family,
                        "organism": organism,
                        "panel_role": panel_role,
                        "source_file": str(path.relative_to(PROJECT_ROOT)),
                        "reviewed_only": row.get("reviewed_only", ""),
                    }
                )
        else:
            for organism in load_species_txt(path):
                key = (family, organism, panel_role)
                if key in seen:
                    continue
                seen.add(key)
                rows.append(
                    {
                        "family": family,
                        "organism": organism,
                        "panel_role": panel_role,
                        "source_file": str(path.relative_to(PROJECT_ROOT)),
                        "reviewed_only": "",
                    }
                )

    rows.sort(key=lambda r: (r["family"], r["panel_role"], r["organism"].lower()))
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build unified species panel metadata table.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = build_panel()
    write_csv(args.output, rows)
    print(f"Wrote {len(rows)} panel entries -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
