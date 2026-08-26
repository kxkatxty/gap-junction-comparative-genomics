#!/usr/bin/env python3
"""Export thesis figures into a single output folder."""

from __future__ import annotations

import argparse
import shutil
import sys
from pathlib import Path

from pipeline.common import PROJECT_ROOT, RESULTS_DIR, SYNVOY_DIR

DEFAULT_OUT = RESULTS_DIR / "figures" / "thesis"


def collect_sources() -> list[tuple[Path, str]]:
    sources: list[tuple[Path, str]] = []
    exon_plots = RESULTS_DIR / "plots"
    if exon_plots.exists():
        for path in sorted(exon_plots.glob("**/*.png")):
            if not path.is_file():
                continue
            rel = path.relative_to(exon_plots)
            sources.append((path, f"exon_{rel.as_posix().replace('/', '_')}"))

    for root in RESULTS_DIR.glob("*_discovery"):
        for path in sorted(root.glob("*/plots/*.png")):
            rel = path.relative_to(root)
            sources.append((path, f"discovery_{root.name}_{rel.as_posix().replace('/', '_')}"))

    synvoy = SYNVOY_DIR / "results" / "innexin_synvoy"
    if synvoy.exists():
        for path in sorted(synvoy.glob("**/*synteny_plot.html")):
            if not path.is_file():
                continue
            rel = path.relative_to(synvoy)
            sources.append((path, f"synvoy_{rel.as_posix().replace('/', '_')}"))

    return sources


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Copy thesis figures to one export folder.")
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT)
    parser.add_argument("--clean", action="store_true", help="Remove output dir before export.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.clean and args.output_dir.exists():
        shutil.rmtree(args.output_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    copied = 0
    for src, name in collect_sources():
        dst = args.output_dir / name
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        copied += 1

    manifest = args.output_dir / "manifest.txt"
    manifest.write_text(
        "\n".join(name for _, name in collect_sources()) + ("\n" if copied else ""),
        encoding="utf-8",
    )
    print(f"Exported {copied} figures -> {args.output_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
