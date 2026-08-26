#!/usr/bin/env python3
"""Summarize candidate validation status across discovery runs."""

from __future__ import annotations

import argparse
import sys
from collections import Counter
from pathlib import Path

from pipeline.common import METADATA_DIR, REJECTED_RANK, write_csv
from pipeline.family_result_merger import collect_candidates

DEFAULT_OUT = METADATA_DIR / "candidate_validation_summary.csv"


def classify_confidence(rank: str) -> str:
    if rank == REJECTED_RANK:
        return "reject"
    if "high_confidence" in rank:
        return "high"
    if "possible" in rank:
        return "possible"
    if "weak" in rank:
        return "weak"
    return "other"


def build_summary() -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    by_species: dict[tuple[str, str, str], list[dict[str, str]]] = {}
    for cand in collect_candidates(accepted_only=False):
        key = (cand["family"], cand["discovery_batch"], cand["organism"])
        by_species.setdefault(key, []).append(cand)

    for (family, batch, organism), cands in sorted(by_species.items()):
        ranks = Counter(classify_confidence(c["rank_category"]) for c in cands)
        rows.append(
            {
                "family": family,
                "discovery_batch": batch,
                "organism": organism,
                "candidate_total": str(len(cands)),
                "accepted_total": str(sum(1 for c in cands if c["rank_category"] != REJECTED_RANK)),
                "high_confidence": str(ranks.get("high", 0)),
                "possible": str(ranks.get("possible", 0)),
                "weak": str(ranks.get("weak", 0)),
                "rejected": str(ranks.get("reject", 0)),
                "best_rank_category": max(
                    (c for c in cands if c["rank_category"] != REJECTED_RANK),
                    key=lambda c: float(c["rank_score"] or 0),
                    default={"rank_category": ""},
                )["rank_category"],
                "best_rank_score": max(
                    (float(c["rank_score"] or 0) for c in cands if c["rank_category"] != REJECTED_RANK),
                    default=0.0,
                ).__str__(),
            }
        )
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize candidate validation by species.")
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    rows = build_summary()
    write_csv(args.output, rows)
    print(f"Wrote validation summary for {len(rows)} species -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
