from __future__ import annotations

import argparse
import csv
import re
import sys
import time
from dataclasses import dataclass, field
from pathlib import Path

from download_gff3_from_species_list import (
    entrez_esearch_assembly,
    entrez_esummary_assembly,
    gff_candidate_urls,
    rank_assemblies,
    url_exists,
)

PROJECT_ROOT = Path("project")
METADATA_DIR = PROJECT_ROOT / "metadata"
DEFAULT_CLADES_CSV = METADATA_DIR / "target_clades_innexin.csv"
DEFAULT_CATALOG_OUT = METADATA_DIR / "species_clade_catalog.csv"
DEFAULT_DISCOVERY_QUEUE = METADATA_DIR / "species_discovery_queue.txt"
DEFAULT_ANNOTATED_QUEUE = METADATA_DIR / "species_annotated_pipeline_queue.txt"
MASTER_CHECKED_CSV = METADATA_DIR / "species_master_checked.csv"
ANNOTATION_SUMMARY_CSV = METADATA_DIR / "family_annotation_summary_new_species.csv"

REQUEST_DELAY_SEC = 0.34
DEFAULT_MODEL_SPECIES: tuple[str, ...] = (
    "Drosophila melanogaster",
    "Caenorhabditis elegans",
    "Caenorhabditis briggsae",
    "Danio rerio",
    "Homo sapiens",
    "Mus musculus",
    "Rattus norvegicus",
    "Gallus gallus",
    "Xenopus laevis",
    "Xenopus tropicalis",
    "Daphnia pulex",
    "Daphnia magna",
    "Hydra vulgaris",
    "Nematostella vectensis",
    "Amphimedon queenslandica",
    "Trichoplax adhaerens",
    "Mnemiopsis leidyi",
)

ASSEMBLY_LEVEL_SCORES = {
    "complete genome": 50,
    "chromosome": 40,
    "scaffold": 25,
    "contig": 10,
}


@dataclass
class CladeTarget:
    clade: str
    ncbi_search_term: str
    notes: str = ""
    priority: int = 99


@dataclass
class SpeciesCandidate:
    clade: str
    organism: str
    assembly_accession: str = ""
    assembly_level: str = ""
    assembly_source: str = ""
    has_gff_annotation: str = "no"
    innexin_annotation_status: str = ""
    connexin_annotation_status: str = ""
    pannexin_annotation_status: str = ""
    already_checked: str = "no"
    innexin_discovery_run: str = "no"
    model_species: str = "no"
    priority_score: int = 0
    priority_tier: str = "low"
    recommended_pipeline: str = "skip"
    clade_priority: int = 99
    notes: str = ""


def normalize_key(name: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", name.lower())


def read_clade_targets(path: Path) -> list[CladeTarget]:
    targets: list[CladeTarget] = []
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            clade = (row.get("clade") or "").strip()
            term = (row.get("ncbi_search_term") or "").strip()
            if not clade or not term:
                continue
            targets.append(
                CladeTarget(
                    clade=clade,
                    ncbi_search_term=term,
                    notes=(row.get("notes") or "").strip(),
                    priority=int(row.get("priority") or 99),
                )
            )
    return sorted(targets, key=lambda item: (item.priority, item.clade.lower()))


def load_master_checked(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {normalize_key(row["organism"]): row for row in rows if row.get("organism")}


def load_annotation_summary(path: Path) -> dict[str, dict[str, str]]:
    if not path.exists():
        return {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    return {normalize_key(row["organism"]): row for row in rows if row.get("organism")}


def assembly_source_label(record: dict[str, str]) -> str:
    if record.get("FtpPath_RefSeq"):
        return "RefSeq"
    if record.get("FtpPath_GenBank"):
        return "GenBank"
    return ""


def assembly_level_score(level: str) -> int:
    lowered = (level or "").lower()
    for key, score in ASSEMBLY_LEVEL_SCORES.items():
        if key in lowered:
            return score
    return 0


def has_gff_for_assembly(record: dict[str, str], *, check_remote: bool = True) -> bool:
    if not check_remote:
        return bool(record.get("FtpPath_RefSeq") or record.get("FtpPath_GenBank"))
    for url, _base, _source in gff_candidate_urls(record):
        if url_exists(url):
            return True
    return False


def clade_assembly_search_term(base_term: str) -> str:
    if "latest[filter]" in base_term.lower():
        return base_term
    return f"{base_term} AND latest[filter]"


def fetch_assemblies_for_clade(
    target: CladeTarget,
    *,
    retmax: int = 500,
    request_delay_sec: float = REQUEST_DELAY_SEC,
) -> list[dict[str, str]]:
    term = clade_assembly_search_term(target.ncbi_search_term)
    uids = entrez_esearch_assembly(target.clade, term, retmax=retmax)
    time.sleep(request_delay_sec)
    if not uids:
        return []
    return entrez_esummary_assembly(uids)


def clean_organism_name(name: str) -> str:
    return re.sub(r"\s*\([^)]+\)\s*$", "", name.strip())


def best_assembly_per_organism(records: list[dict[str, str]]) -> dict[str, dict[str, str]]:
    grouped: dict[str, list[dict[str, str]]] = {}
    for record in records:
        organism = clean_organism_name(
            (record.get("Organism") or record.get("SpeciesName") or "").strip()
        )
        if not organism:
            continue
        grouped.setdefault(organism, []).append(record)

    best: dict[str, dict[str, str]] = {}
    for organism, items in grouped.items():
        ranked = rank_assemblies(items)
        if ranked:
            best[organism] = ranked[0]
    return best


def innexin_is_annotated(status: str) -> bool:
    return status == "annotated_family_member_found"


def family_is_annotated(status: str) -> bool:
    return status == "annotated_family_member_found"


def compute_priority(candidate: SpeciesCandidate, *, family: str = "innexin") -> SpeciesCandidate:
    score = 0
    annotation_status = {
        "innexin": candidate.innexin_annotation_status,
        "connexin": candidate.connexin_annotation_status,
        "pannexin": candidate.pannexin_annotation_status,
    }.get(family, candidate.innexin_annotation_status)

    if candidate.assembly_accession:
        score += 20
    score += assembly_level_score(candidate.assembly_level)
    if candidate.has_gff_annotation == "yes":
        score += 15
    if candidate.model_species == "yes":
        score -= 25

    if annotation_status == "annotated_family_member_found":
        candidate.recommended_pipeline = "annotated"
        score -= 20
    elif annotation_status == "related_or_family_like_entry_found":
        candidate.recommended_pipeline = "review_related"
        score += 5
    elif annotation_status == "no_related_entry_found":
        candidate.recommended_pipeline = "discovery"
        score += 30
    elif candidate.already_checked == "no":
        candidate.recommended_pipeline = "classify_then_discovery"
        score += 25
    else:
        candidate.recommended_pipeline = "classify_first"
        score += 15

    if candidate.innexin_discovery_run == "yes":
        score -= 10

    candidate.priority_score = score
    if score >= 55:
        candidate.priority_tier = "high"
    elif score >= 35:
        candidate.priority_tier = "medium"
    else:
        candidate.priority_tier = "low"
    return candidate


def build_candidate(
    target: CladeTarget,
    organism: str,
    assembly: dict[str, str],
    *,
    master: dict[str, dict[str, str]],
    summary: dict[str, dict[str, str]],
    model_species: set[str],
    check_gff: bool,
    family: str = "innexin",
) -> SpeciesCandidate:
    key = normalize_key(organism)
    master_row = master.get(key, {})
    summary_row = summary.get(key, {})

    candidate = SpeciesCandidate(
        clade=target.clade,
        organism=organism,
        assembly_accession=assembly.get("AssemblyAccession", ""),
        assembly_level=assembly.get("AssemblyStatus", ""),
        assembly_source=assembly_source_label(assembly),
        has_gff_annotation="yes" if has_gff_for_assembly(assembly, check_remote=check_gff) else "no",
        innexin_annotation_status=summary_row.get("innexin_category", ""),
        connexin_annotation_status=summary_row.get("connexin_category", ""),
        pannexin_annotation_status=summary_row.get("pannexin_category", ""),
        already_checked=master_row.get("annotation_searched", "no") or "no",
        innexin_discovery_run=master_row.get("innexin_discovery_run", "no") or "no",
        model_species="yes" if key in model_species else "no",
        clade_priority=target.priority,
        notes=target.notes,
    )
    return compute_priority(candidate, family=family)


def collect_candidates_for_clade(
    target: CladeTarget,
    *,
    master: dict[str, dict[str, str]],
    summary: dict[str, dict[str, str]],
    model_species: set[str],
    check_gff: bool,
    retmax: int,
    request_delay_sec: float,
    family: str = "innexin",
) -> list[SpeciesCandidate]:
    records = fetch_assemblies_for_clade(
        target,
        retmax=retmax,
        request_delay_sec=request_delay_sec,
    )
    best_by_organism = best_assembly_per_organism(records)
    candidates = [
        build_candidate(
            target,
            organism,
            assembly,
            master=master,
            summary=summary,
            model_species=model_species,
            check_gff=check_gff,
            family=family,
        )
        for organism, assembly in sorted(best_by_organism.items(), key=lambda item: item[0].lower())
    ]
    candidates.sort(key=lambda item: (-item.priority_score, item.organism.lower()))
    return candidates


def filter_candidates(
    candidates: list[SpeciesCandidate],
    *,
    exclude_model_species: bool,
    exclude_already_checked: bool,
    min_assembly_level: str | None,
    require_gff: bool,
    pipeline: str | None,
) -> list[SpeciesCandidate]:
    min_score = assembly_level_score(min_assembly_level or "")
    kept: list[SpeciesCandidate] = []
    for candidate in candidates:
        if exclude_model_species and candidate.model_species == "yes":
            continue
        if exclude_already_checked and candidate.already_checked == "yes":
            continue
        if require_gff and candidate.has_gff_annotation != "yes":
            continue
        if min_assembly_level and assembly_level_score(candidate.assembly_level) < min_score:
            continue
        if pipeline and candidate.recommended_pipeline != pipeline:
            continue
        kept.append(candidate)
    return kept


def apply_per_clade_limit(
    candidates: list[SpeciesCandidate],
    per_clade_limit: int | None,
) -> list[SpeciesCandidate]:
    if per_clade_limit is None:
        return candidates
    grouped: dict[str, list[SpeciesCandidate]] = {}
    for candidate in candidates:
        grouped.setdefault(candidate.clade, []).append(candidate)
    limited: list[SpeciesCandidate] = []
    for clade in sorted(grouped):
        items = sorted(
            grouped[clade],
            key=lambda item: (-item.priority_score, item.organism.lower()),
        )
        limited.extend(items[:per_clade_limit])
    return limited


def write_catalog(path: Path, candidates: list[SpeciesCandidate]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "clade",
        "organism",
        "assembly_accession",
        "assembly_level",
        "assembly_source",
        "has_gff_annotation",
        "innexin_annotation_status",
        "connexin_annotation_status",
        "pannexin_annotation_status",
        "already_checked",
        "innexin_discovery_run",
        "model_species",
        "priority_score",
        "priority_tier",
        "recommended_pipeline",
        "clade_priority",
        "notes",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for candidate in candidates:
            writer.writerow(
                {
                    "clade": candidate.clade,
                    "organism": candidate.organism,
                    "assembly_accession": candidate.assembly_accession,
                    "assembly_level": candidate.assembly_level,
                    "assembly_source": candidate.assembly_source,
                    "has_gff_annotation": candidate.has_gff_annotation,
                    "innexin_annotation_status": candidate.innexin_annotation_status,
                    "connexin_annotation_status": candidate.connexin_annotation_status,
                    "pannexin_annotation_status": candidate.pannexin_annotation_status,
                    "already_checked": candidate.already_checked,
                    "innexin_discovery_run": candidate.innexin_discovery_run,
                    "model_species": candidate.model_species,
                    "priority_score": candidate.priority_score,
                    "priority_tier": candidate.priority_tier,
                    "recommended_pipeline": candidate.recommended_pipeline,
                    "clade_priority": candidate.clade_priority,
                    "notes": candidate.notes,
                }
            )


def write_species_list(path: Path, organisms: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    unique = list(dict.fromkeys(organisms))
    path.write_text("\n".join(unique) + ("\n" if unique else ""), encoding="utf-8")


def summarize_by_clade(candidates: list[SpeciesCandidate]) -> list[dict[str, str]]:
    grouped: dict[str, list[SpeciesCandidate]] = {}
    for candidate in candidates:
        grouped.setdefault(candidate.clade, []).append(candidate)

    rows: list[dict[str, str]] = []
    for clade in sorted(grouped):
        items = grouped[clade]
        rows.append(
            {
                "clade": clade,
                "species_with_assembly": str(len(items)),
                "annotated_innexin": str(
                    sum(1 for item in items if innexin_is_annotated(item.innexin_annotation_status))
                ),
                "no_innexin_annotation": str(
                    sum(
                        1
                        for item in items
                        if item.innexin_annotation_status == "no_related_entry_found"
                    )
                ),
                "discovery_recommended": str(
                    sum(1 for item in items if item.recommended_pipeline == "discovery")
                ),
                "high_priority": str(sum(1 for item in items if item.priority_tier == "high")),
            }
        )
    return rows


def write_clade_summary(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "clade",
        "species_with_assembly",
        "annotated_innexin",
        "no_innexin_annotation",
        "discovery_recommended",
        "high_priority",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Collect species systematically by clade from NCBI Assembly, "
            "filter by assembly/annotation quality, and rank candidates for "
            "annotation or innexin discovery pipelines."
        )
    )
    parser.add_argument(
        "--family",
        choices=("innexin", "connexin", "pannexin"),
        default="innexin",
        help="Gene family used for pipeline priority and discovery queue (default: innexin).",
    )
    parser.add_argument(
        "--clades-csv",
        type=Path,
        default=DEFAULT_CLADES_CSV,
        help=f"CSV with target clades (default: {DEFAULT_CLADES_CSV})",
    )
    parser.add_argument(
        "--clade",
        action="append",
        default=[],
        help="Run only these clade names (repeatable). Default: all clades in CSV.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=DEFAULT_CATALOG_OUT,
        help=f"Species catalog CSV (default: {DEFAULT_CATALOG_OUT})",
    )
    parser.add_argument(
        "--clade-summary",
        type=Path,
        default=METADATA_DIR / "species_clade_summary.csv",
        help="Per-clade summary CSV.",
    )
    parser.add_argument(
        "--discovery-queue",
        type=Path,
        default=DEFAULT_DISCOVERY_QUEUE,
        help="Species list for innexin discovery pipeline.",
    )
    parser.add_argument(
        "--annotated-queue",
        type=Path,
        default=DEFAULT_ANNOTATED_QUEUE,
        help="Species list for annotated-protein pipeline.",
    )
    parser.add_argument(
        "--per-clade-limit",
        type=int,
        default=10,
        help="Keep only the top N species per clade after ranking (default: 10).",
    )
    parser.add_argument(
        "--retmax",
        type=int,
        default=500,
        help="Maximum NCBI assembly records to fetch per clade (default: 500).",
    )
    parser.add_argument(
        "--exclude-model-species",
        action="store_true",
        help="Drop common model species from the output.",
    )
    parser.add_argument(
        "--exclude-already-checked",
        action="store_true",
        help="Drop species already present in species_master_checked.csv.",
    )
    parser.add_argument(
        "--min-assembly-level",
        choices=["contig", "scaffold", "chromosome", "complete genome"],
        default="scaffold",
        help="Minimum assembly level to keep (default: scaffold).",
    )
    parser.add_argument(
        "--require-gff",
        action="store_true",
        help="Keep only species whose best assembly has a reachable GFF/GFF3 file.",
    )
    parser.add_argument(
        "--skip-gff-check",
        action="store_true",
        help="Do not probe NCBI FTP for GFF files (faster; marks has_gff_annotation=no).",
    )
    parser.add_argument(
        "--include-model-species-list",
        action="append",
        default=[],
        help="Additional model species to exclude when --exclude-model-species is set.",
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

    master = load_master_checked(MASTER_CHECKED_CSV)
    summary = load_annotation_summary(ANNOTATION_SUMMARY_CSV)
    model_species = {normalize_key(name) for name in DEFAULT_MODEL_SPECIES}
    model_species.update(normalize_key(name) for name in args.include_model_species_list)

    all_candidates: list[SpeciesCandidate] = []
    for index, target in enumerate(targets, start=1):
        print(f"[{index}/{len(targets)}] Collecting assemblies for {target.clade}")
        candidates = collect_candidates_for_clade(
            target,
            master=master,
            summary=summary,
            model_species=model_species,
            check_gff=not args.skip_gff_check,
            retmax=args.retmax,
            request_delay_sec=REQUEST_DELAY_SEC,
            family=args.family,
        )
        print(f"  -> {len(candidates)} species with assemblies")
        all_candidates.extend(candidates)

    filtered = filter_candidates(
        all_candidates,
        exclude_model_species=args.exclude_model_species,
        exclude_already_checked=args.exclude_already_checked,
        min_assembly_level=args.min_assembly_level,
        require_gff=args.require_gff,
        pipeline=None,
    )
    filtered.sort(key=lambda item: (item.clade_priority, -item.priority_score, item.organism.lower()))
    filtered = apply_per_clade_limit(filtered, args.per_clade_limit)

    write_catalog(args.output, filtered)
    write_clade_summary(args.clade_summary, summarize_by_clade(filtered))

    discovery_queue = list(
        dict.fromkeys(
            item.organism
            for item in filtered
            if item.recommended_pipeline
            in {"discovery", "classify_then_discovery", "classify_first"}
        )
    )
    annotated_queue = [item.organism for item in filtered if item.recommended_pipeline == "annotated"]
    write_species_list(args.discovery_queue, discovery_queue)
    write_species_list(args.annotated_queue, annotated_queue)

    print(f"Saved catalog: {args.output} ({len(filtered)} species)")
    print(f"Saved clade summary: {args.clade_summary}")
    print(f"Saved discovery queue: {args.discovery_queue} ({len(discovery_queue)} species)")
    print(f"Saved annotated queue: {args.annotated_queue} ({len(annotated_queue)} species)")
    return 0


if __name__ == "__main__":
    sys.exit(main())
