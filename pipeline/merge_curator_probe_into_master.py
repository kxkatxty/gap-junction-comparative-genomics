#!/usr/bin/env python3
"""Merge gene-curator present/fragmentary loci into master candidate + panel tables."""

from __future__ import annotations

import argparse
import csv
from pathlib import Path

import pandas as pd

from pipeline.common import METADATA_DIR, RESULTS_DIR

CURATOR_SUMMARY = RESULTS_DIR / "gene_curator_probe" / "innexin_search_summary.csv"
CURATOR_ROOT = RESULTS_DIR / "gene_curator_probe"
CAND_PATH = METADATA_DIR / "gap_junction_candidates_master.csv"
PANEL_PATH = METADATA_DIR / "species_panel_master.csv"
BATCH = "gene_curator_probe"


def _organism(species_dir: str) -> str:
    return species_dir.replace("_", " ")


def _load_hits(species: str) -> list[dict[str, str]]:
    path = CURATOR_ROOT / species / "hits.csv"
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def curator_rows(*, include_known: bool) -> list[dict[str, object]]:
    summary = list(csv.DictReader(CURATOR_SUMMARY.open(encoding="utf-8")))
    rows: list[dict[str, object]] = []
    idx = 0
    for rec in summary:
        if rec.get("status") != "searched":
            continue
        species = rec["species"]
        prior = rec.get("prior_discovery") or ""
        try:
            prior_accepted = int(rec["prior_accepted"]) if rec.get("prior_accepted") not in ("", None) else None
        except ValueError:
            prior_accepted = None
        is_new = prior == "not_run" or prior_accepted == 0
        if not include_known and not is_new:
            continue
        for hit in _load_hits(species):
            if hit.get("verdict") not in {"present", "fragmentary"}:
                continue
            idx += 1
            length = int(hit["aa_length"])
            verdict = hit["verdict"]
            # Unify with curator classifier present gate (280 aa). Do not invent
            # rubric-style scores: leave rank_score empty and mark evidence tier.
            if verdict == "present" and length >= 280:
                rank = "curator_present_innexin"
            elif verdict == "present":
                rank = "weak_manual_review"
            else:
                rank = "weak_manual_review"
            rows.append(
                {
                    "family": "innexin",
                    "discovery_batch": BATCH,
                    "species_dir": species,
                    "organism": _organism(species),
                    "candidate_id": f"curator_inx_{idx:03d}",
                    "seqid": hit["seqid"],
                    "strand": hit["strand"],
                    "locus_start": int(hit["start"]),
                    "locus_end": int(hit["end"]),
                    "exon_count": "",
                    "protein_length": length,
                    "tm_helix_count": "",
                    "reference_identity": float(hit["identity"]),
                    "reference_coverage": "",
                    "best_reference_hit": hit.get("target", ""),
                    "rank_category": rank,
                    "rank_score": "",  # not a discovery rubric score
                    "validation_flags": f"curator_{verdict};evidence_tier=threshold_classifier",
                    "results_dir": f"project/results/gene_curator_probe/{species}",
                }
            )
    return rows


def merge_candidates(*, include_known: bool) -> tuple[int, int, int]:
    cand = pd.read_csv(CAND_PATH)
    before = len(cand)
    cand = cand[cand["discovery_batch"] != BATCH].copy()
    new_rows = curator_rows(include_known=include_known)
    if new_rows:
        cand = pd.concat([cand, pd.DataFrame(new_rows)], ignore_index=True)
    cand.to_csv(CAND_PATH, index=False)
    accepted = int((cand["rank_category"] != "rejected_false_positive").sum())
    return before, len(cand), accepted


def merge_panel() -> int:
    panel = pd.read_csv(PANEL_PATH)
    existing = {str(o).lower().replace(" ", "_").rstrip(".") for o in panel["organism"]}
    summary = list(csv.DictReader(CURATOR_SUMMARY.open(encoding="utf-8")))
    added = 0
    for rec in summary:
        if int(rec.get("present_loci") or 0) <= 0:
            continue
        species = rec["species"]
        key = species.lower().rstrip(".")
        if key in existing or any(key == e or key in e or e in key for e in existing):
            # still add exact missing Adineta_vaga / Rotaria_macrura
            org = _organism(species)
            if any(str(o).lower() == org.lower() for o in panel["organism"]):
                continue
        panel = pd.concat(
            [
                panel,
                pd.DataFrame(
                    [
                        {
                            "family": "innexin",
                            "organism": _organism(species),
                            "panel_role": "discovery_innexin_curator",
                            "source_file": "project/results/gene_curator_probe/innexin_search_summary.csv",
                            "reviewed_only": "",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
        existing.add(key)
        added += 1
    panel.to_csv(PANEL_PATH, index=False)
    return added


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--include-known",
        action="store_true",
        help="Also merge curator loci for species that already had accepted discovery hits",
    )
    return p.parse_args()


def main() -> int:
    args = parse_args()
    before, after, accepted = merge_candidates(include_known=args.include_known)
    added_panel = merge_panel()
    print(f"Candidates: {before} → {after} (accepted non-rejected: {accepted})")
    print(f"Panel species added: {added_panel} (total {len(pd.read_csv(PANEL_PATH))})")
    print(f"Wrote {CAND_PATH}")
    print(f"Wrote {PANEL_PATH}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
