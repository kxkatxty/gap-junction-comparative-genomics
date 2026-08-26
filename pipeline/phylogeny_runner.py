#!/usr/bin/env python3
"""Run IQ-TREE on a protein alignment."""

from __future__ import annotations

import argparse
import shutil
import subprocess
import sys
from pathlib import Path

from pipeline.common import RESULTS_DIR

DEFAULT_OUT_DIR = RESULTS_DIR / "phylogeny"


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build phylogenetic tree with IQ-TREE.")
    parser.add_argument("alignment", type=Path, help="Aligned FASTA (MAFFT output).")
    parser.add_argument(
        "--prefix",
        type=Path,
        default=None,
        help="Output prefix (default: alignment stem in results/phylogeny/).",
    )
    parser.add_argument(
        "--model",
        default="MFP",
        help="IQ-TREE model (-m), default MFP (ModelFinder Plus).",
    )
    parser.add_argument("--threads", type=int, default=2)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.alignment.exists():
        print(f"Alignment not found: {args.alignment}", file=sys.stderr)
        return 1
    if shutil.which("iqtree") is None and shutil.which("iqtree2") is None:
        print("iqtree/iqtree2 not found on PATH", file=sys.stderr)
        return 1

    prefix = args.prefix or (DEFAULT_OUT_DIR / args.alignment.stem)
    prefix.parent.mkdir(parents=True, exist_ok=True)
    exe = "iqtree2" if shutil.which("iqtree2") else "iqtree"

    cmd = [
        exe,
        "-s",
        str(args.alignment),
        "-m",
        args.model,
        "-bb",
        "1000",
        "-nt",
        str(args.threads),
        "-pre",
        str(prefix),
    ]
    subprocess.run(cmd, check=True)
    tree = Path(f"{prefix}.treefile")
    print(f"Tree written -> {tree}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
