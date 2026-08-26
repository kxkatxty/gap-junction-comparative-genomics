#!/usr/bin/env python3
"""Merge discovery candidate tables across species into one master dataset."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pipeline.common import (
    DISCOVERY_DIRS,
    METADATA_DIR,
    REJECTED_RANK,
    infer_family_from_path,
    read_csv_rows,
    write_csv,
)

DEFAULT_OUT = METADATA_DIR / "gap_junction_candidates_master.csv"


def collect_candidates(*, accepted_only: bool = False) -> list[dict[str, str]]:
    merged: list[dict[str, str]] = []
    for root in DISCOVERY_DIRS:
        if not root.exists():
            continue
        family = infer_family_from_path(root)
        batch = root.name
        for csv_path in sorted(root.glob("*/candidate_loci.csv")):
            organism_dir = csv_path.parent.name
            for row in read_csv_rows(csv_path):
                if accepted_only and row.get("rank_category") == REJECTED_RANK:
                    continue
                merged.append(
                    {
                        "family": family,
                        "discovery_batch": batch,
                        "species_dir": organism_dir,
                        "organism": row.get("organism", organism_dir.replace("_", " ")),
                        "candidate_id": row.get("candidate_id", ""),
                        "seqid": row.get("seqid", ""),
                        "strand": row.get("strand", ""),
                        "locus_start": row.get("locus_start", ""),
                        "locus_end": row.get("locus_end", ""),
                        "exon_count": row.get("exon_count", ""),
                        "protein_length": row.get("protein_length", ""),
                        "tm_helix_count": row.get("tm_helix_count", ""),
                        "reference_identity": row.get("reference_identity", ""),
                        "reference_coverage": row.get("reference_coverage", ""),
                        "best_reference_hit": row.get("best_reference_hit", ""),
                        "rank_category": row.get("rank_category", ""),
                        "rank_score": row.get("rank_score", ""),
                        "validation_flags": row.get("validation_flags", ""),
                        "results_dir": str(csv_path.parent.relative_to(root.parent.parent)),
                    }
                )
    return merged


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Merge discovery candidate tables.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    parser.add_argument(
        "--accepted-only",
        action="store_true",
        help="Keep only non-rejected candidates.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = collect_candidates(accepted_only=args.accepted_only)
    write_csv(args.output, rows)
    species = len({r["organism"] for r in rows})
    print(f"Wrote {len(rows)} candidates from {species} species -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
