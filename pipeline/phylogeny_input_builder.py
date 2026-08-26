#!/usr/bin/env python3
"""Build a multi-FASTA for phylogeny from references and/or accepted discovery candidates."""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

from pipeline.common import REFERENCES_DIR, RESULTS_DIR, REJECTED_RANK
from pipeline.family_result_merger import collect_candidates

DEFAULT_OUT = RESULTS_DIR / "phylogeny" / "gap_junction_phylogeny_input.fasta"


def read_fasta_records(path: Path) -> list[tuple[str, str]]:
    if not path.exists():
        return []
    records: list[tuple[str, str]] = []
    header = ""
    seq_parts: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        if line.startswith(">"):
            if header:
                records.append((header, "".join(seq_parts)))
            header = line[1:].strip()
            seq_parts = []
        else:
            seq_parts.append(line.strip())
    if header:
        records.append((header, "".join(seq_parts)))
    return records


def write_fasta(path: Path, records: list[tuple[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines: list[str] = []
    for header, seq in records:
        lines.append(f">{header}")
        for i in range(0, len(seq), 80):
            lines.append(seq[i : i + 80])
    path.write_text("\n".join(lines) + ("\n" if lines else ""), encoding="utf-8")


def build_records(
    *,
    family: str | None,
    include_references: bool,
    accepted_only: bool,
) -> list[tuple[str, str]]:
    records: list[tuple[str, str]] = []
    seen: set[str] = set()

    if include_references:
        ref_root = REFERENCES_DIR / f"{family}s" if family else REFERENCES_DIR
        for fasta in sorted(ref_root.rglob("*.fasta")):
            for header, seq in read_fasta_records(fasta):
                key = header.split()[0]
                if key in seen:
                    continue
                seen.add(key)
                records.append((header, seq))

    for cand in collect_candidates(accepted_only=accepted_only):
        if family and cand["family"] != family:
            continue
        if accepted_only and cand["rank_category"] == REJECTED_RANK:
            continue
        fasta_path = RESULTS_DIR / cand["discovery_batch"] / cand["species_dir"] / (
            "candidate_innexins.fasta" if cand["family"] == "innexin" else "candidate_connexins.fasta"
        )
        if not fasta_path.exists():
            continue
        for header, seq in read_fasta_records(fasta_path):
            if cand["candidate_id"] not in header and cand["candidate_id"] not in header.replace(" ", "_"):
                continue
            label = f"{cand['organism']}|{cand['candidate_id']}|{cand['family']}"
            if label in seen:
                continue
            seen.add(label)
            records.append((label, seq))
    return records


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build phylogeny input FASTA.")
    parser.add_argument("--family", choices=("innexin", "connexin"), default=None)
    parser.add_argument("--include-references", action="store_true")
    parser.add_argument("--accepted-only", action="store_true", default=True)
    parser.add_argument("--output", type=Path, default=DEFAULT_OUT)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    records = build_records(
        family=args.family,
        include_references=args.include_references,
        accepted_only=args.accepted_only,
    )
    if not records:
        print("No sequences found for phylogeny input.", file=sys.stderr)
        return 1
    write_fasta(args.output, records)
    print(f"Wrote {len(records)} sequences -> {args.output}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
