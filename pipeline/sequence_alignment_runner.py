#!/usr/bin/env python3
"""Run MAFFT alignment on a multi-FASTA protein file."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from pipeline.common import RESULTS_DIR

DEFAULT_OUT_DIR = RESULTS_DIR / "phylogeny"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Align protein sequences with MAFFT.")
    parser.add_argument("input_fasta", type=Path, help="Multi-sequence FASTA input.")
    parser.add_argument(
        "--output",
        type=Path,
        default=None,
        help="Aligned FASTA output (default: <input>.aln.fasta).",
    )
    parser.add_argument(
        "--algorithm",
        choices=("auto", "linsi", "ginsi", "fftns"),
        default="auto",
        help="MAFFT mode (auto uses --auto).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.input_fasta.exists():
        print(f"Input not found: {args.input_fasta}", file=sys.stderr)
        return 1
    if shutil.which("mafft") is None:
        print("mafft not found on PATH", file=sys.stderr)
        return 1

    out = args.output or args.input_fasta.with_suffix(".aln.fasta")
    out.parent.mkdir(parents=True, exist_ok=True)

    if args.algorithm == "auto":
        mafft_args = ["mafft", "--auto", str(args.input_fasta)]
    else:
        mafft_args = ["mafft", f"--{args.algorithm}", str(args.input_fasta)]

    result = subprocess.run(
        mafft_args,
        check=True,
        capture_output=True,
        text=True,
    )
    out.write_text(result.stdout, encoding="utf-8")
    print(f"Wrote alignment -> {out}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
