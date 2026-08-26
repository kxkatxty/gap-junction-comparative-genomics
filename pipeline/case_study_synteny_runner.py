#!/usr/bin/env python3
"""Run SynVoy for one or more selected innexin/connexin case-study genes."""

from __future__ import annotations

import argparse
import csv
import subprocess
import sys
from pathlib import Path

from pipeline.common import METADATA_DIR, PROJECT_ROOT, SYNVOY_DIR

DEFAULT_QUEUE = METADATA_DIR / "innexin_synvoy_queue.tsv"
LOG_DIR = PROJECT_ROOT / "project" / "results" / "synvoy_innexins" / "logs"

SYNVOY_FAST_ARGS = [
    "--n_flanking_genes",
    "5",
    "--adaptive_max_regions",
    "3",
    "--enable_smith_waterman",
    "false",
    "--exon_level_search",
    "false",
]


def load_queue(path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            rows.append(row)
    return rows


def job_done(outdir_name: str) -> bool:
    outdir = SYNVOY_DIR / "results" / "innexin_synvoy" / outdir_name
    return any(outdir.glob("synteny_block_*_synteny_plot.html"))


def run_job(row: dict[str, str], *, max_genomes: int) -> int:
    accession = row["accession"]
    gene = row["gene"]
    outdir_name = row["outdir_name"]
    outdir = f"results/innexin_synvoy/{outdir_name}"
    log_path = LOG_DIR / f"{outdir_name}.log"
    LOG_DIR.mkdir(parents=True, exist_ok=True)

    if job_done(outdir_name):
        print(f"SKIP {gene} ({accession}) already complete")
        return 0

    cmd = [
        "nextflow",
        "run",
        "main.nf",
        "-profile",
        "standard",
        "--mode",
        "easy",
        "--query_id",
        accession,
        "--max_genomes",
        str(max_genomes),
        "--outdir",
        outdir,
        "--auto_params",
        "false",
        "--multi_profile",
        "false",
        *SYNVOY_FAST_ARGS,
        "-resume",
    ]
    print(f"START {gene} ({accession}) -> {outdir}")
    with log_path.open("w", encoding="utf-8") as log_handle:
        proc = subprocess.run(cmd, cwd=SYNVOY_DIR, stdout=log_handle, stderr=subprocess.STDOUT)
    if proc.returncode == 0:
        print(f"OK {gene} ({accession})")
    else:
        print(f"FAIL {gene} ({accession}) exit={proc.returncode} (see {log_path})")
    return proc.returncode


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Run SynVoy case study for selected genes.")
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--gene", action="append", default=[], help="Gene symbol(s) to run.")
    parser.add_argument("--accession", action="append", default=[], help="UniProt accession(s).")
    parser.add_argument("--max-genomes", type=int, default=3)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    queue = load_queue(args.queue)
    if not queue:
        print(f"No jobs in queue: {args.queue}", file=sys.stderr)
        return 1

    selected = queue
    if args.gene or args.accession:
        wanted_genes = {g.casefold() for g in args.gene}
        wanted_acc = {a.upper() for a in args.accession}
        selected = [
            row
            for row in queue
            if row["gene"].casefold() in wanted_genes or row["accession"].upper() in wanted_acc
        ]

    if not selected:
        print("No matching jobs in queue.", file=sys.stderr)
        return 1

    failures = 0
    for row in selected:
        failures += run_job(row, max_genomes=args.max_genomes) != 0
    return 1 if failures else 0


if __name__ == "__main__":
    sys.exit(main())
