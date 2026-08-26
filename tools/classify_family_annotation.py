from __future__ import annotations

import argparse
import csv
import json
import re
import sys
import time
import xml.etree.ElementTree as ET
from dataclasses import dataclass, field
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


USER_AGENT = "BachelorCursor-family-classifier/1.0"
PROJECT_ROOT = Path("project")
METADATA_DIR = PROJECT_ROOT / "metadata"
DEFAULT_OUTPUT = METADATA_DIR / "family_annotation_status.csv"
DEFAULT_SUMMARY = METADATA_DIR / "family_annotation_summary.csv"

FAMILIES = ("innexin", "connexin", "pannexin")

CATEGORY_ANNOTATED = "annotated_family_member_found"
CATEGORY_RELATED = "related_or_family_like_entry_found"
CATEGORY_NONE = "no_related_entry_found"

CONFIDENCE_HIGH = "high"
CONFIDENCE_MEDIUM = "medium"

REQUEST_DELAY_SEC = 0.34

FAMILY_SYMBOL_PREFIXES: dict[str, tuple[str, ...]] = {
    "connexin": ("GJA", "GJB", "GJC", "GJD", "GJE", "GJF"),
    "pannexin": ("PANX",),
    "innexin": ("INX",),
}

INNEXIN_MODEL_ALIASES: tuple[str, ...] = (
    "unc-7",
    "unc-9",
    "eat-5",
    "shakb",
    "shak-b",
    "ogre",
    "zpg",
)

LEVEL_B_TERMS: dict[str, tuple[str, ...]] = {
    "connexin": (
        "gap junction",
        "gap junction protein",
        "connexin-like",
        "junctional channel",
    ),
    "pannexin": (
        "gap junction",
        "gap junction protein",
        "pannexin-like",
        "innexin homolog",
    ),
    "innexin": (
        "gap junction",
        "gap junction protein",
        "innexin-like",
        "innexin homolog",
    ),
}

LOW_CONFIDENCE_PHRASES: tuple[str, ...] = (
    "channel protein",
    "membrane protein",
    "junctional membrane protein",
    "transmembrane protein",
    "ion channel",
)


@dataclass
class SearchHit:
    database: str
    level: str
    query: str
    accession: str = ""
    gene_name: str = ""
    protein_name: str = ""
    description: str = ""
    extra: dict[str, str] = field(default_factory=dict)


@dataclass
class Classification:
    final_category: str
    confidence: str
    matching_term: str
    database_source: str
    level: str
    hit: SearchHit | None = None


def load_species_from_csv(csv_path: Path) -> list[str]:
    organisms: list[str] = []
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        for row in reader:
            organism = (row.get("organism") or "").strip()
            if organism:
                organisms.append(organism)
    return dedupe_preserve_order(organisms)


def load_species_from_text(path: Path) -> list[str]:
    names: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            names.append(line)
    return dedupe_preserve_order(names)


def dedupe_preserve_order(items: list[str]) -> list[str]:
    seen: set[str] = set()
    out: list[str] = []
    for item in items:
        key = item.casefold()
        if key not in seen:
            seen.add(key)
            out.append(item)
    return out


def http_get_json(url: str, timeout: int = 60) -> dict:
    request = Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    with urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def http_get_text(url: str, timeout: int = 60) -> str:
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=timeout) as response:
        return response.read().decode("utf-8")


def organism_clause_uniprot(organism: str) -> str:
    return f'organism_name:"{organism}"'


def organism_clause_ncbi(organism: str) -> str:
    return f'"{organism}"[Organism]'


def normalize_token(text: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (text or "").lower())


def gene_matches_prefix(gene_name: str, prefix: str) -> bool:
    token = normalize_token(gene_name)
    prefix_norm = normalize_token(prefix)
    if not token or not prefix_norm:
        return False
    return token.startswith(prefix_norm)


def gene_matches_alias(gene_name: str, alias: str) -> bool:
    return normalize_token(gene_name) == normalize_token(alias)


def combined_text(hit: SearchHit) -> str:
    return " ".join(
        part
        for part in (hit.protein_name, hit.gene_name, hit.description, hit.accession)
        if part
    ).lower()


def contains_family_word(text: str, family: str) -> bool:
    pattern = rf"\b{re.escape(family.lower())}\b"
    return bool(re.search(pattern, text.lower()))


def contains_phrase(text: str, phrase: str) -> bool:
    return phrase.lower() in text.lower()


def has_characteristic_symbol(family: str, hit: SearchHit) -> str | None:
    gene = hit.gene_name or ""
    for prefix in FAMILY_SYMBOL_PREFIXES.get(family, ()):
        if gene_matches_prefix(gene, prefix):
            return gene or prefix
    if family == "innexin":
        for alias in INNEXIN_MODEL_ALIASES:
            if gene_matches_alias(gene, alias):
                return gene or alias
    return None


def is_low_confidence_only(text: str, family: str) -> bool:
    lowered = text.lower()
    if contains_family_word(lowered, family):
        return False
    if any(term in lowered for term in ("gap junction", "innexin-like", "connexin-like", "pannexin-like")):
        return False
    if any(contains_phrase(lowered, phrase) for phrase in LOW_CONFIDENCE_PHRASES):
        return True
    return False


def classify_hit(family: str, hit: SearchHit) -> Classification | None:
    text = combined_text(hit)

    if contains_family_word(text, family):
        return Classification(
            final_category=CATEGORY_ANNOTATED,
            confidence=CONFIDENCE_HIGH,
            matching_term=family,
            database_source=hit.database,
            level=hit.level,
            hit=hit,
        )

    symbol_match = has_characteristic_symbol(family, hit)
    if symbol_match:
        return Classification(
            final_category=CATEGORY_ANNOTATED,
            confidence=CONFIDENCE_HIGH,
            matching_term=symbol_match,
            database_source=hit.database,
            level=hit.level,
            hit=hit,
        )

    if hit.level == "B":
        for term in LEVEL_B_TERMS.get(family, ()):
            if contains_phrase(text, term):
                return Classification(
                    final_category=CATEGORY_RELATED,
                    confidence=CONFIDENCE_MEDIUM,
                    matching_term=term,
                    database_source=hit.database,
                    level=hit.level,
                    hit=hit,
                )

        if family == "innexin":
            for alias in INNEXIN_MODEL_ALIASES:
                if contains_phrase(text, alias) or gene_matches_alias(hit.gene_name, alias):
                    return Classification(
                        final_category=CATEGORY_RELATED,
                        confidence=CONFIDENCE_MEDIUM,
                        matching_term=alias,
                        database_source=hit.database,
                        level=hit.level,
                        hit=hit,
                    )

        if family == "pannexin" and contains_phrase(text, "innexin homolog"):
            return Classification(
                final_category=CATEGORY_RELATED,
                confidence=CONFIDENCE_MEDIUM,
                matching_term="innexin homolog",
                database_source=hit.database,
                level=hit.level,
                hit=hit,
            )

    if is_low_confidence_only(text, family):
        return None

    if hit.level == "B" and text.strip():
        return Classification(
            final_category=CATEGORY_RELATED,
            confidence=CONFIDENCE_MEDIUM,
            matching_term=hit.protein_name or hit.gene_name or hit.description,
            database_source=hit.database,
            level=hit.level,
            hit=hit,
        )

    return None


def pick_best_classification(candidates: list[Classification]) -> Classification | None:
    if not candidates:
        return None

    priority = {
        CATEGORY_ANNOTATED: 0,
        CATEGORY_RELATED: 1,
        CATEGORY_NONE: 2,
    }
    confidence_rank = {CONFIDENCE_HIGH: 0, CONFIDENCE_MEDIUM: 1, "": 2}
    level_rank = {"A": 0, "B": 1}

    return sorted(
        candidates,
        key=lambda item: (
            priority[item.final_category],
            confidence_rank.get(item.confidence, 9),
            level_rank.get(item.level, 9),
            item.matching_term,
        ),
    )[0]


def merge_family_results(
  level_a: list[Classification],
  level_b: list[Classification],
) -> Classification:
    best = pick_best_classification(level_a + level_b)
    if best:
        return best
    return Classification(
        final_category=CATEGORY_NONE,
        confidence="",
        matching_term="",
        database_source="",
        level="",
        hit=None,
    )


def uniprot_search(query: str, size: int = 5) -> tuple[int, list[dict]]:
    url = (
        "https://rest.uniprot.org/uniprotkb/search?"
        f"query={quote(query)}&format=json&size={size}"
        "&fields=accession,gene_names,protein_name,organism_name,cc_function,ft_domain,keyword"
    )
    try:
        payload = http_get_json(url)
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
        return 0, []
    total = int(payload.get("total", 0) or 0)
    return total, payload.get("results", [])


def uniprot_result_to_hit(
    entry: dict, *, database: str, level: str, query: str
) -> SearchHit:
    genes = entry.get("genes", []) or []
    gene_name = ""
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")

    protein_desc = entry.get("proteinDescription", {})
    protein_name = ""
    rec = protein_desc.get("recommendedName", {})
    if rec:
        protein_name = rec.get("fullName", {}).get("value", "")
    if not protein_name:
        subs = protein_desc.get("submissionNames", []) or []
        if subs:
            protein_name = subs[0].get("fullName", {}).get("value", "")

    cc_bits: list[str] = []
    for cc in entry.get("comments", []) or []:
        if cc.get("commentType") == "FUNCTION":
            texts = cc.get("texts", []) or []
            for text in texts:
                value = text.get("value", "")
                if value:
                    cc_bits.append(value)

    keywords = ", ".join(
        item.get("name", "")
        for item in entry.get("keywords", []) or []
        if item.get("name")
    )

    description = "; ".join(part for part in (", ".join(cc_bits), keywords) if part)

    return SearchHit(
        database=database,
        level=level,
        query=query,
        accession=entry.get("primaryAccession", ""),
        gene_name=gene_name,
        protein_name=protein_name,
        description=description,
    )


def ncbi_esearch(db: str, term: str, retmax: int = 5) -> tuple[int, list[str]]:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db={db}&term={quote(term)}&retmode=json&retmax={retmax}"
    )
    try:
        payload = http_get_json(url)
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
        return 0, []
    result = payload.get("esearchresult", {})
    count = int(result.get("count", 0) or 0)
    ids = result.get("idlist", []) or []
    return count, ids


def ncbi_gene_summary(gene_id: str) -> dict[str, str]:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=gene&id={gene_id}&retmode=json"
    )
    try:
        payload = http_get_json(url)
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
        return {}
    result = payload.get("result", {})
    if gene_id not in result:
        return {}
    rec = result[gene_id]
    return {
        "accession": gene_id,
        "gene_name": rec.get("name", "") or rec.get("nomenclaturesymbol", ""),
        "protein_name": rec.get("description", "") or rec.get("summary", ""),
        "description": rec.get("summary", "") or rec.get("description", ""),
    }


def ncbi_protein_summary(protein_id: str) -> dict[str, str]:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=protein&id={protein_id}&retmode=json"
    )
    try:
        payload = http_get_json(url)
    except (HTTPError, URLError, TimeoutError, json.JSONDecodeError):
        return {}
    result = payload.get("result", {})
    if protein_id not in result:
        return {}
    rec = result[protein_id]
    return {
        "accession": rec.get("accessionversion", "") or protein_id,
        "gene_name": rec.get("title", "").split()[0] if rec.get("title") else "",
        "protein_name": rec.get("title", ""),
        "description": rec.get("title", ""),
    }


def ncbi_ids_to_hits(
    db: str,
    ids: list[str],
    *,
    level: str,
    query: str,
) -> list[SearchHit]:
    hits: list[SearchHit] = []
    for item_id in ids:
        if db == "gene":
            summary = ncbi_gene_summary(item_id)
        else:
            summary = ncbi_protein_summary(item_id)
        if not summary:
            continue
        hits.append(
            SearchHit(
                database=f"ncbi_{db}",
                level=level,
                query=query,
                accession=summary.get("accession", item_id),
                gene_name=summary.get("gene_name", ""),
                protein_name=summary.get("protein_name", ""),
                description=summary.get("description", ""),
            )
        )
        time.sleep(REQUEST_DELAY_SEC)
    return hits


def level_a_queries(family: str, organism: str) -> list[tuple[str, str]]:
    org_u = organism_clause_uniprot(organism)
    org_n = organism_clause_ncbi(organism)
    queries: list[tuple[str, str]] = []

    queries.append(
        (
            "uniprot",
            f"{org_u} AND protein_name:{family}",
        )
    )
    queries.append(
        (
            "uniprot",
            f'{org_u} AND (gene:{family} OR annotation:{family})',
        )
    )
    for prefix in FAMILY_SYMBOL_PREFIXES.get(family, ()):
        queries.append(("uniprot", f"{org_u} AND gene:{prefix}*"))

    queries.append(("ncbi_gene", f"{family}[All Fields] AND {org_n}"))
    queries.append(("ncbi_gene", f"{family}[Gene Name] AND {org_n}"))

    for prefix in FAMILY_SYMBOL_PREFIXES.get(family, ()):
        queries.append(("ncbi_gene", f"{prefix}*[Gene Symbol] AND {org_n}"))

    if family == "innexin":
        for alias in INNEXIN_MODEL_ALIASES:
            queries.append(("ncbi_gene", f"{alias}[Gene Symbol] AND {org_n}"))

    return queries


def level_b_queries(family: str, organism: str) -> list[tuple[str, str]]:
    org_u = organism_clause_uniprot(organism)
    org_n = organism_clause_ncbi(organism)
    queries: list[tuple[str, str]] = []

    for term in LEVEL_B_TERMS.get(family, ()):
        if " " in term:
            queries.append(("uniprot", f'{org_u} AND (protein_name:"{term}" OR cc_function:"{term}")'))
            queries.append(("ncbi_gene", f'"{term}"[All Fields] AND {org_n}'))
        else:
            queries.append(("uniprot", f"{org_u} AND (protein_name:{term} OR cc_function:{term})"))
            queries.append(("ncbi_gene", f"{term}[All Fields] AND {org_n}"))

    for prefix in FAMILY_SYMBOL_PREFIXES.get(family, ()):
        if ("uniprot", f"{org_u} AND gene:{prefix}*") not in queries:
            queries.append(("uniprot", f"{org_u} AND gene:{prefix}*"))
        queries.append(("ncbi_gene", f"{prefix}*[Gene Symbol] AND {org_n}"))

    if family == "innexin":
        for alias in INNEXIN_MODEL_ALIASES:
            queries.append(("uniprot", f'{org_u} AND gene_exact:{alias}'))
            queries.append(("ncbi_gene", f"{alias}[Gene Symbol] AND {org_n}"))

    return dedupe_query_pairs(queries)


def dedupe_query_pairs(queries: list[tuple[str, str]]) -> list[tuple[str, str]]:
    seen: set[tuple[str, str]] = set()
    out: list[tuple[str, str]] = []
    for db, query in queries:
        key = (db, query)
        if key in seen:
            continue
        seen.add(key)
        out.append(key)
    return out


def run_query(db: str, query: str, *, level: str) -> list[SearchHit]:
    if db == "uniprot":
        _total, results = uniprot_search(query)
        hits = [
            uniprot_result_to_hit(entry, database="uniprot", level=level, query=query)
            for entry in results
        ]
        time.sleep(REQUEST_DELAY_SEC)
        return hits

    if db == "ncbi_gene":
        _count, ids = ncbi_esearch("gene", query)
        time.sleep(REQUEST_DELAY_SEC)
        return ncbi_ids_to_hits("gene", ids, level=level, query=query)

    if db == "ncbi_protein":
        _count, ids = ncbi_esearch("protein", query)
        time.sleep(REQUEST_DELAY_SEC)
        return ncbi_ids_to_hits("protein", ids, level=level, query=query)

    return []


def run_level(family: str, organism: str, level: str) -> tuple[list[Classification], int]:
    if level == "A":
        query_pairs = level_a_queries(family, organism)
    else:
        query_pairs = level_b_queries(family, organism)

    classifications: list[Classification] = []
    raw_hit_count = 0

    for db, query in query_pairs:
        hits = run_query(db, query, level=level)
        raw_hit_count += len(hits)
        for hit in hits:
            result = classify_hit(family, hit)
            if result:
                classifications.append(result)

    return classifications, raw_hit_count


def run_protein_fallback(family: str, organism: str) -> list[Classification]:
    org_n = organism_clause_ncbi(organism)
    queries = [
        (f"{family}[Title] AND {org_n}", "A"),
        (f'"gap junction"[Title] AND {org_n}', "B"),
    ]
    for term in LEVEL_B_TERMS.get(family, ())[:2]:
        queries.append((f'"{term}"[Title] AND {org_n}', "B"))

    out: list[Classification] = []
    for query, level in queries:
        hits = run_query("ncbi_protein", query, level=level)
        for hit in hits:
            result = classify_hit(family, hit)
            if result:
                out.append(result)
    return out


def classify_species_family(organism: str, family: str) -> dict[str, str]:
    level_a, level_a_count = run_level(family, organism, "A")
    best_a = pick_best_classification(level_a)

    level_b: list[Classification] = []
    level_b_count = 0
    if not best_a or best_a.final_category != CATEGORY_ANNOTATED:
        level_b, level_b_count = run_level(family, organism, "B")
        if not pick_best_classification(level_b):
            level_b.extend(run_protein_fallback(family, organism))

    final = merge_family_results(level_a, level_b)
    hit = final.hit

    return {
        "organism": organism,
        "family": family,
        "exact_annotation_found": "yes"
        if final.final_category == CATEGORY_ANNOTATED
        else "no",
        "related_annotation_found": "yes"
        if final.final_category == CATEGORY_RELATED
        else "no",
        "final_category": final.final_category,
        "confidence": final.confidence,
        "top_matching_term": final.matching_term,
        "database_source": final.database_source,
        "search_level": final.level,
        "top_accession": hit.accession if hit else "",
        "top_gene_symbol": hit.gene_name if hit else "",
        "top_protein_name": hit.protein_name if hit else "",
        "level_a_hit_count": str(level_a_count),
        "level_b_hit_count": str(level_b_count),
        "example_query": hit.query if hit else "",
    }


def write_results(rows: list[dict[str, str]], out_path: Path) -> None:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "organism",
        "family",
        "exact_annotation_found",
        "related_annotation_found",
        "final_category",
        "confidence",
        "top_matching_term",
        "database_source",
        "search_level",
        "top_accession",
        "top_gene_symbol",
        "top_protein_name",
        "level_a_hit_count",
        "level_b_hit_count",
        "example_query",
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


def write_summary(rows: list[dict[str, str]], out_path: Path) -> None:
    by_organism: dict[str, dict[str, str]] = {}
    for row in rows:
        organism = row["organism"]
        by_organism.setdefault(organism, {"organism": organism})
        family = row["family"]
        by_organism[organism][f"{family}_category"] = row["final_category"]
        by_organism[organism][f"{family}_term"] = row["top_matching_term"]

    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["organism"] + [
        f"{family}_{suffix}"
        for family in FAMILIES
        for suffix in ("category", "term")
    ]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for organism in sorted(by_organism):
            writer.writerow(by_organism[organism])


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Classify species x family annotation status using layered UniProt "
            "and NCBI searches."
        )
    )
    parser.add_argument(
        "--species-csv",
        type=Path,
        help="CSV with an organism column (e.g. species_config.csv).",
    )
    parser.add_argument(
        "-f",
        "--file",
        action="append",
        default=[],
        type=Path,
        help="Text file with one species per line.",
    )
    parser.add_argument(
        "species",
        nargs="*",
        help="Species names on the command line.",
    )
    parser.add_argument(
        "-o",
        "--output",
        type=Path,
        default=DEFAULT_OUTPUT,
        help=f"Detailed CSV output (default: {DEFAULT_OUTPUT}).",
    )
    parser.add_argument(
        "--summary",
        type=Path,
        default=DEFAULT_SUMMARY,
        help=f"Wide summary CSV (default: {DEFAULT_SUMMARY}).",
    )
    parser.add_argument(
        "--families",
        nargs="+",
        choices=FAMILIES,
        default=list(FAMILIES),
        help="Families to search (default: all three).",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N species (useful for testing).",
    )
    return parser.parse_args(argv)


def resolve_species(args: argparse.Namespace) -> list[str]:
    organisms: list[str] = []
    if args.species_csv:
        organisms.extend(load_species_from_csv(args.species_csv))
    for path in args.file:
        organisms.extend(load_species_from_text(path))
    organisms.extend(args.species)
    organisms = dedupe_preserve_order(organisms)
    if args.limit is not None:
        organisms = organisms[: args.limit]
    return organisms


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    organisms = resolve_species(args)
    if not organisms:
        print("No species provided. Use --species-csv, --file, or positional species names.")
        return 1

    rows: list[dict[str, str]] = []
    total = len(organisms) * len(args.families)
    step = 0

    print(f"Classifying {len(organisms)} species x {len(args.families)} families...")
    for organism in organisms:
        for family in args.families:
            step += 1
            print(f"[{step}/{total}] {organism} / {family}")
            row = classify_species_family(organism, family)
            rows.append(row)
            print(
                f"  -> {row['final_category']}"
                + (f" ({row['top_matching_term']})" if row["top_matching_term"] else "")
            )

    write_results(rows, args.output)
    write_summary(rows, args.summary)
    print(f"Saved detailed results to {args.output}")
    print(f"Saved summary to {args.summary}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
