#!/usr/bin/env python3
"""Collect species by evolutionary clade until N without annotated connexin are found."""

from __future__ import annotations

import argparse
import csv
import sys
import time
from pathlib import Path

from classify_family_annotation import CATEGORY_ANNOTATED, classify_species_family
from collect_species_by_clade import (
    DEFAULT_MODEL_SPECIES,
    CladeTarget,
    REQUEST_DELAY_SEC,
    assembly_level_score,
    assembly_source_label,
    best_assembly_per_organism,
    clean_organism_name,
    fetch_assemblies_for_clade,
    normalize_key,
    read_clade_targets,
)

PROJECT_ROOT = Path("project")
METADATA_DIR = PROJECT_ROOT / "metadata"
DEFAULT_CLADES = METADATA_DIR / "target_clades_connexin_evolutionary.csv"
DEFAULT_CATALOG = METADATA_DIR / "species_clade_catalog_connexin_evolutionary.csv"
DEFAULT_QUEUE = METADATA_DIR / "not_annotated_connexin_evolutionary.txt"
DEFAULT_STATUS = METADATA_DIR / "family_annotation_status_connexin_evolutionary.csv"


def connexin_not_annotated(category: str) -> bool:
    return category != CATEGORY_ANNOTATED


def load_exclude_organisms(paths: list[Path]) -> set[str]:
    excluded: set[str] = set()
    for path in paths:
        if not path.exists():
            continue
        if path.suffix == ".csv":
            with path.open("r", encoding="utf-8", newline="") as handle:
                for row in csv.DictReader(handle):
                    organism = row.get("organism", "").strip()
                    if organism:
                        excluded.add(normalize_key(organism))
            continue
        for line in path.read_text(encoding="utf-8").splitlines():
            line = line.strip()
            if line and not line.startswith("#"):
                excluded.add(normalize_key(line))
    return excluded


def rank_organisms(
    best_by_organism: dict[str, dict[str, str]],
    *,
    model_species: set[str],
    min_assembly_level: str,
) -> list[tuple[str, dict[str, str], int]]:
    min_score = assembly_level_score(min_assembly_level)
    ranked: list[tuple[str, dict[str, str], int]] = []
    for organism, assembly in best_by_organism.items():
        if normalize_key(organism) in model_species:
            continue
        level = assembly.get("AssemblyStatus", "")
        if assembly_level_score(level) < min_score:
            continue
        score = 20 + assembly_level_score(level)
        ranked.append((organism, assembly, score))
    ranked.sort(key=lambda item: (-item[2], item[0].lower()))
    return ranked


def collect_clade_targets(
    target: CladeTarget,
    *,
    per_clade_limit: int,
    model_species: set[str],
    exclude_organisms: set[str],
    min_assembly_level: str,
    retmax: int,
    request_delay_sec: float,
) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    records = fetch_assemblies_for_clade(
        target, retmax=retmax, request_delay_sec=request_delay_sec
    )
    ranked_organisms = rank_organisms(
        best_assembly_per_organism(records),
        model_species=model_species,
        min_assembly_level=min_assembly_level,
    )

    selected: list[dict[str, str]] = []
    screened: list[dict[str, str]] = []

    for organism, assembly, assembly_score in ranked_organisms:
        if len(selected) >= per_clade_limit:
            break
        if normalize_key(organism) in exclude_organisms:
            print(f"  skip {organism} (already used)", flush=True)
            continue

        print(f"  screening {organism} ...", flush=True)
        row = classify_species_family(organism, "connexin")
        category = row["final_category"]
        screened.append(row)
        print(f"    -> {category}" + (f" ({row['top_matching_term']})" if row["top_matching_term"] else ""))

        if not connexin_not_annotated(category):
            continue

        selected.append(
            {
                "clade": target.clade,
                "organism": organism,
                "assembly_accession": assembly.get("AssemblyAccession", ""),
                "assembly_level": assembly.get("AssemblyStatus", ""),
                "assembly_source": assembly_source_label(assembly),
                "connexin_category": category,
                "connexin_term": row.get("top_matching_term", ""),
                "connexin_evidence_level": row.get("final_evidence_level", ""),
                "connexin_source": row.get("final_source", ""),
                "assembly_priority_score": str(assembly_score),
                "notes": target.notes,
            }
        )
        time.sleep(request_delay_sec)

    return selected, screened


def write_catalog(path: Path, rows: list[dict[str, str]]) -> None:
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    fieldnames = list(rows[0].keys())
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_species_list(path: Path, organisms: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(organisms) + ("\n" if organisms else ""), encoding="utf-8")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Sample species by evolutionary clade until N species without "
            "annotated connexin are found per clade."
        )
    )
    parser.add_argument("--clades-csv", type=Path, default=DEFAULT_CLADES)
    parser.add_argument("--output", type=Path, default=DEFAULT_CATALOG)
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--status", type=Path, default=DEFAULT_STATUS)
    parser.add_argument(
        "--per-clade-limit",
        type=int,
        default=5,
        help="Stop after this many species without annotated connexin per clade.",
    )
    parser.add_argument("--retmax", type=int, default=500)
    parser.add_argument(
        "--min-assembly-level",
        choices=["contig", "scaffold", "chromosome", "complete genome"],
        default="scaffold",
    )
    parser.add_argument("--clade", action="append", default=[])
    parser.add_argument(
        "--exclude-txt",
        type=Path,
        action="append",
        default=[],
        help="Species list(s) to skip (e.g. prior batch queue).",
    )
    parser.add_argument(
        "--exclude-catalog",
        type=Path,
        action="append",
        default=[],
        help="Prior catalog CSV(s) whose organisms should be skipped.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    targets = read_clade_targets(args.clades_csv)
    if args.clade:
        wanted = {name.casefold() for name in args.clade}
        targets = [target for target in targets if target.clade.casefold() in wanted]
    if not targets:
        print("No clade targets found.")
        return 1

    model_species = {normalize_key(name) for name in DEFAULT_MODEL_SPECIES}
    exclude_organisms = load_exclude_organisms(args.exclude_txt + args.exclude_catalog)
    all_selected: list[dict[str, str]] = []
    all_screened: list[dict[str, str]] = []

    for index, target in enumerate(targets, start=1):
        print(
            f"[{index}/{len(targets)}] {target.clade}: "
            f"collect until {args.per_clade_limit} without annotated connexin"
            + (f" (excluding {len(exclude_organisms)} prior species)" if exclude_organisms else "")
        )
        selected, screened = collect_clade_targets(
            target,
            per_clade_limit=args.per_clade_limit,
            model_species=model_species,
            exclude_organisms=exclude_organisms,
            min_assembly_level=args.min_assembly_level,
            retmax=args.retmax,
            request_delay_sec=REQUEST_DELAY_SEC,
        )
        all_selected.extend(selected)
        all_screened.extend(screened)
        print(f"  -> selected {len(selected)}/{args.per_clade_limit} (screened {len(screened)})")

    write_catalog(args.output, all_selected)
    write_species_list(args.queue, [row["organism"] for row in all_selected])
    if all_screened:
        fieldnames = list(all_screened[0].keys())
        with args.status.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=fieldnames)
            writer.writeheader()
            writer.writerows(all_screened)

    print(f"Saved catalog: {args.output} ({len(all_selected)} species)")
    print(f"Saved queue: {args.queue}")
    print(f"Saved screened status: {args.status} ({len(all_screened)} rows)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
