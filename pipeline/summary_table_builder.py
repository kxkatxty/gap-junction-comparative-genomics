#!/usr/bin/env python3
"""Build the master thesis summary table (panel + annotation + discovery + quality)."""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

from pipeline.candidate_validator import build_summary as build_validation_summary
from pipeline.common import METADATA_DIR, PROJECT_ROOT, read_csv_rows, write_csv
from pipeline.species_panel_builder import build_panel

DEFAULT_OUT = METADATA_DIR / "thesis_master_summary.csv"


def load_annotation_lookup() -> dict[tuple[str, str], dict[str, str]]:
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    for path in sorted(METADATA_DIR.glob("family_annotation_status*.csv")):
        for row in read_csv_rows(path):
            organism = row.get("organism", "").strip()
            family = row.get("family", "").strip().lower()
            if not organism or not family:
                continue
            key = (family, organism)
            if key not in lookup:
                lookup[key] = row
    return lookup


def load_quality_lookup() -> dict[str, dict[str, str]]:
    path = METADATA_DIR / "gap_junction_quality_table.csv"
    lookup: dict[str, dict[str, str]] = {}
    for row in read_csv_rows(path):
        organism = row.get("organism", row.get("species", "")).strip()
        if not organism:
            continue
        if organism not in lookup:
            lookup[organism] = row
    return lookup


def load_discovery_lookup() -> dict[tuple[str, str], dict[str, str]]:
    lookup: dict[tuple[str, str], dict[str, str]] = {}
    for row in build_validation_summary():
        lookup[(row["family"], row["organism"])] = row
    return lookup


def build_master_summary() -> list[dict[str, str]]:
    annotation = load_annotation_lookup()
    quality = load_quality_lookup()
    discovery = load_discovery_lookup()
    rows: list[dict[str, str]] = []

    for panel in build_panel():
        family = panel["family"]
        organism = panel["organism"]
        ann = annotation.get((family, organism), {})
        qual = quality.get(organism, {})
        disc = discovery.get((family, organism), {})

        rows.append(
            {
                "family": family,
                "organism": organism,
                "panel_role": panel["panel_role"],
                "panel_source": panel["source_file"],
                "annotation_category": ann.get("final_category", ann.get("connexin_category", "")),
                "annotation_term": ann.get("top_matching_term", ann.get("connexin_term", "")),
                "annotation_evidence": ann.get("final_evidence_level", ""),
                "assembly_level": qual.get("assembly_level", qual.get("assembly_status", "")),
                "has_gff_annotation": qual.get("annotation_available", qual.get("has_gff_annotation", "")),
                "discovery_candidates": disc.get("candidate_total", "0"),
                "discovery_accepted": disc.get("accepted_total", "0"),
                "discovery_high_confidence": disc.get("high_confidence", "0"),
                "discovery_best_rank": disc.get("best_rank_category", ""),
                "discovery_best_score": disc.get("best_rank_score", ""),
            }
        )
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build master thesis summary table.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = build_master_summary()
    write_csv(args.output, rows)
    print(f"Wrote {len(rows)} master summary rows -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
