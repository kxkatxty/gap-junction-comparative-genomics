from __future__ import annotations

import argparse
import csv
import json
import sys
import time
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen

USER_AGENT = "gap-junction-config-builder/1.0"
SUPPORTED_FAMILIES = ("innexin", "connexin", "pannexin")
DEFAULT_FAMILIES = ("innexin", "connexin", "pannexin")
VERTEBRATA_TAXON_ID = 7742


def load_species_names(paths: list[Path], inline: list[str]) -> list[str]:
    names: list[str] = []
    for item in inline:
        for part in item.replace(";", ",").split(","):
            name = part.strip()
            if name:
                names.append(name)

    for path in paths:
        text = path.read_text(encoding="utf-8")
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            names.append(line)

    seen: set[str] = set()
    unique: list[str] = []
    for name in names:
        key = name.casefold()
        if key not in seen:
            seen.add(key)
            unique.append(name)
    return unique


def uniprot_get(url: str) -> dict:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def search_uniprot_count(query: str) -> int:
    url = (
        "https://rest.uniprot.org/uniprotkb/search?"
        f"query={quote(query)}&format=json&size=0"
    )
    payload = uniprot_get(url)
    return int(payload.get("total", 0) or 0)


def build_protein_query(family: str, organism: str, reviewed_only: bool) -> str:
    terms = [f"protein_name:{family}", f'organism_name:"{organism}"']
    if reviewed_only:
        terms.append("reviewed:true")
    return " AND ".join(terms)


def _lineage_has_vertebrata(lineage: list) -> bool:
    for node in lineage:
        if not isinstance(node, dict):
            continue
        if node.get("taxonId") == VERTEBRATA_TAXON_ID:
            return True
        if node.get("scientificName") == "Vertebrata":
            return True
    return False


def _pick_taxonomy_match(organism: str, results: list) -> dict | None:
    target = organism.casefold().strip()
    for entry in results:
        if entry.get("scientificName", "").casefold() == target:
            return entry
    for entry in results:
        synonyms = entry.get("synonyms", []) or []
        if any(str(s).casefold() == target for s in synonyms):
            return entry
    return results[0] if results else None


def taxonomy_is_vertebrate(organism: str) -> bool | None:
    queries = [
        f'scientific_name:"{organism}"',
        f"scientific_name:{organism}",
        organism,
    ]
    for taxonomy_query in queries:
        url = (
            "https://rest.uniprot.org/taxonomy/search?"
            f"query={quote(taxonomy_query)}&format=json&size=10"
        )
        try:
            payload = uniprot_get(url)
        except (HTTPError, URLError):
            continue

        results = payload.get("results", [])
        match = _pick_taxonomy_match(organism, results)
        if not match:
            continue

        lineage = match.get("lineage", [])
        if _lineage_has_vertebrata(lineage):
            return True

        # Metazoa without Vertebrata → invertebrate (innexin family).
        is_metazoa = any(
            isinstance(node, dict) and node.get("scientificName") == "Metazoa"
            for node in lineage
        )
        if is_metazoa:
            return False

    return None


def detect_family(
    organism: str,
    reviewed_only: bool,
    families_to_try: tuple[str, ...],
) -> tuple[str, str]:
    """
    Return (family, method) where method explains how the family was chosen.
    """
    is_vertebrate = taxonomy_is_vertebrate(organism)
    time.sleep(0.1)
    if is_vertebrate is True:
        return "connexin", "taxonomy"
    if is_vertebrate is False:
        return "innexin", "taxonomy"

    hits: list[tuple[str, int]] = []
    for family in families_to_try:
        query = build_protein_query(family, organism, reviewed_only)
        count = search_uniprot_count(query)
        if count > 0:
            hits.append((family, count))
        time.sleep(0.1)

    if len(hits) == 1:
        return hits[0][0], "uniprot"

    if len(hits) > 1:
        priority = {family: index for index, family in enumerate(SUPPORTED_FAMILIES)}
        hits.sort(key=lambda item: (priority[item[0]], -item[1]))
        return hits[0][0], "uniprot_multiple"

    return "connexin", "default_fallback"


def write_config_csv(
    rows: list[dict[str, str]],
    output_path: Path,
) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    with output_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["family", "organism", "reviewed_only"],
        )
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Build a download config CSV (family,organism,reviewed_only) "
            "from a list of species names."
        ),
    )
    parser.add_argument(
        "species",
        nargs="*",
        help="Species names (e.g. 'Homo sapiens' 'Danio rerio').",
    )
    parser.add_argument(
        "-f",
        "--file",
        action="append",
        default=[],
        type=Path,
        help="Text file with one species per line (# comments allowed).",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=Path("project/metadata/species_config.csv"),
        help="Output CSV path (default: project/metadata/species_config.csv).",
    )
    parser.add_argument(
        "--family",
        choices=SUPPORTED_FAMILIES,
        help="Force one family for every species (skip auto-detection).",
    )
    parser.add_argument(
        "--reviewed-only",
        choices=("true", "false"),
        default="true",
        help="Value for reviewed_only column (default: true).",
    )
    parser.add_argument(
        "--no-auto-detect",
        action="store_true",
        help="Require --family; do not query UniProt for family assignment.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    species_names = load_species_names(args.file, args.species)

    if not species_names:
        print(
            "No species provided. Pass names on the command line and/or use --file.",
            file=sys.stderr,
        )
        return 1

    if args.no_auto_detect and not args.family:
        print("--no-auto-detect requires --family.", file=sys.stderr)
        return 1

    reviewed_value = args.reviewed_only
    reviewed_only = reviewed_value == "true"
    rows: list[dict[str, str]] = []

    for organism in species_names:
        if args.family:
            family = args.family
            method = "forced"
        elif args.no_auto_detect:
            family = args.family or "connexin"
            method = "forced"
        else:
            family, method = detect_family(
                organism,
                reviewed_only=reviewed_only,
                families_to_try=DEFAULT_FAMILIES,
            )

        rows.append(
            {
                "family": family,
                "organism": organism,
                "reviewed_only": reviewed_value,
            }
        )
        print(f"{organism} -> {family} ({method})")

    write_config_csv(rows, args.output)
    print(f"\nWrote {len(rows)} row(s) to {args.output}")
    print("Use with: python download_by_family.py", args.output)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
