from __future__ import annotations

import csv
import gzip
import re
import sys
from collections import defaultdict
from pathlib import Path
from statistics import mean

import pandas as pd


PROJECT_ROOT = Path("project")
REFERENCES_DIR = PROJECT_ROOT / "data" / "references"
ANNOTATION_DIR = PROJECT_ROOT / "data" / "annotations"
RESULTS_DIR = PROJECT_ROOT / "results" / "exon_structures"

RESULTS_DIR.mkdir(parents=True, exist_ok=True)

TRANSCRIPT_FEATURES = {"mrna", "transcript"}
EXCLUDED_GENE_KEYWORDS = {"cnst", "ric1"}  # connexin-associated, not connexin genes
EXCLUDED_GENE_PATTERNS = re.compile(
    r"readthrough|lncrna|pseudogene|mycbp",
    re.IGNORECASE,
)
REFERENCE_GENE_ALIASES: dict[str, list[str]] = {
    "zpg": ["inx4", "zpg"],
    "shakb": ["shaking-b", "shaking b", "shakb", "shak-b"],
}


def normalize_species_name(name: str) -> str:
    return name.strip().replace(" ", "_")


def normalize_gene_symbol(symbol: str) -> str:
    return re.sub(r"[^a-z0-9]", "", (symbol or "").lower())


def family_to_ref_folder(family: str) -> Path:
    family_norm = family.lower()
    if family_norm in {"innexin", "innexins"}:
        return REFERENCES_DIR / "innexins"
    if family_norm in {"connexin", "connexins"}:
        return REFERENCES_DIR / "connexins"
    return REFERENCES_DIR / "pannexins"


def parse_attributes(attr_str: str) -> dict[str, str]:
    attr_str = attr_str.strip()
    attrs: dict[str, str] = {}

    if "=" in attr_str:
        for field in attr_str.split(";"):
            field = field.strip()
            if not field or "=" not in field:
                continue
            key, value = field.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
    else:
        for match in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
            attrs[match.group(1)] = match.group(2)

    return attrs


def safe_int(value: str) -> int:
    return int(float(value))


def load_species_csv(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"family", "organism", "reviewed_only"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        return [{key: (val or "").strip() for key, val in row.items()} for row in reader]


def read_fasta_header(path: Path) -> str:
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                return line.strip()
    return ""


def parse_reference_fasta(path: Path) -> dict[str, str]:
    header = read_fasta_header(path)
    match = re.match(r"^(.+?)__([A-Z0-9]{6,10})\.fasta$", path.name, re.IGNORECASE)
    gene_symbol = ""
    accession = ""
    if match:
        gene_symbol = match.group(1)
        accession = match.group(2).upper()
        if normalize_gene_symbol(gene_symbol) == normalize_gene_symbol(accession):
            gene_symbol = ""

    gn_match = re.search(r"\bGN=(\S+)", header)
    if gn_match:
        gene_symbol = gene_symbol or gn_match.group(1)

    protein_name = ""
    if header.startswith(">"):
        parts = header[1:].split(" OS=")
        if parts:
            protein_name = parts[0].split("|", 2)[-1].strip()

    return {
        "fasta_path": str(path),
        "fasta_header": header,
        "gene_symbol": gene_symbol,
        "uniprot_accession": accession,
        "protein_name": protein_name,
    }


def load_reference_proteins(family: str, organism: str) -> list[dict[str, str]]:
    folder = family_to_ref_folder(family) / normalize_species_name(organism)
    if not folder.exists():
        return []
    return [parse_reference_fasta(path) for path in sorted(folder.glob("*.fasta"))]


def candidate_gene_symbols(reference: dict[str, str]) -> list[str]:
    symbols: list[str] = []
    for value in (
        reference.get("gene_symbol", ""),
        reference.get("uniprot_accession", ""),
    ):
        if value:
            symbols.append(value)
            for alias in REFERENCE_GENE_ALIASES.get(normalize_gene_symbol(value), []):
                symbols.append(alias)

    header = reference.get("fasta_header", "")
    gn_match = re.search(r"\bGN=(\S+)", header)
    if gn_match:
        symbols.append(gn_match.group(1))

    uniprot_id_match = re.search(r"\|([A-Z0-9]+)_", header)
    if uniprot_id_match:
        token = uniprot_id_match.group(1)
        symbols.append(token)
        if token.startswith(("CXA", "CXB", "CXC", "CXD")):
            # CXD (delta) connexins are often annotated as GJC in genome GFF.
            mapping = {"CXA": "GJA", "CXB": "GJB", "CXC": "GJC", "CXD": "GJC"}
            suffix = re.sub(r"[^0-9]", "", token[3:])
            if suffix:
                symbols.append(mapping[token[:3]] + suffix)

    seen: set[str] = set()
    unique: list[str] = []
    for symbol in symbols:
        norm = normalize_gene_symbol(symbol)
        if norm and norm not in seen:
            seen.add(norm)
            unique.append(symbol)
    return unique


def find_annotation_file(family: str, organism: str) -> Path | None:
    base = ANNOTATION_DIR / family.lower() / normalize_species_name(organism)
    if not base.exists():
        return None

    candidates: list[Path] = []
    for pattern in ("*.gff", "*.gff3", "*.gtf", "*.gff.gz", "*.gff3.gz"):
        candidates.extend(base.glob(pattern))

    if not candidates:
        return None

    def sort_key(path: Path) -> tuple[int, str]:
        name = path.name.lower()
        if name.endswith(".gff3.gz"):
            priority = 0
        elif name.endswith(".gff3"):
            priority = 1
        elif name.endswith(".gff.gz"):
            priority = 2
        elif name.endswith(".gff"):
            priority = 3
        else:
            priority = 4
        return priority, path.name

    return sorted(candidates, key=sort_key)[0]


def open_annotation(path: Path):
    if path.suffix.lower() == ".gz" or path.name.lower().endswith(".gff.gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def parse_annotation(annotation_path: Path):
    genes: dict[str, dict] = {}
    transcripts: dict[str, dict] = {}
    exons_by_transcript: dict[str, list[tuple[int, int]]] = defaultdict(list)
    gene_to_transcripts: dict[str, list[str]] = defaultdict(list)

    print(f"  Parsing annotation: {annotation_path.name}")
    with open_annotation(annotation_path) as handle:
        for line_no, line in enumerate(handle, start=1):
            if line_no % 2_000_000 == 0:
                print(f"    ...{line_no:,} lines")

            if not line.strip() or line.startswith("#"):
                continue

            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue

            seqid, _source, feature_type, start, end, _score, strand, _phase, attributes = parts
            attrs = parse_attributes(attributes)
            start_i = safe_int(start)
            end_i = safe_int(end)
            feature = feature_type.lower()

            if feature == "gene":
                gene_id = attrs.get("ID") or attrs.get("gene_id")
                if not gene_id:
                    continue

                gene_symbol = (
                    attrs.get("gene_name")
                    or attrs.get("Name")
                    or attrs.get("gene")
                    or attrs.get("locus_tag")
                    or gene_id
                )
                gene_name = (
                    attrs.get("product")
                    or attrs.get("description")
                    or attrs.get("Note")
                    or gene_symbol
                )
                gene_synonym = attrs.get("gene_synonym", "")
                gene_biotype = attrs.get("gene_biotype", "")

                genes[gene_id] = {
                    "gene_id": gene_id,
                    "gene_symbol": gene_symbol,
                    "gene_name": gene_name,
                    "gene_synonym": gene_synonym,
                    "gene_biotype": gene_biotype,
                    "seqid": seqid,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                }

            elif feature in TRANSCRIPT_FEATURES:
                transcript_id = attrs.get("ID") or attrs.get("transcript_id")
                parent_gene = attrs.get("Parent") or attrs.get("gene_id")
                if not transcript_id:
                    continue

                transcript_name = (
                    attrs.get("Name") or attrs.get("transcript_name") or transcript_id
                )
                transcript_product = attrs.get("product", "")
                transcripts[transcript_id] = {
                    "transcript_id": transcript_id,
                    "parent_gene": parent_gene,
                    "transcript_name": transcript_name,
                    "transcript_product": transcript_product,
                    "seqid": seqid,
                    "start": start_i,
                    "end": end_i,
                    "strand": strand,
                }
                if parent_gene:
                    gene_to_transcripts[parent_gene].append(transcript_id)

            elif feature == "exon":
                parent = attrs.get("Parent") or attrs.get("transcript_id")
                if not parent:
                    continue
                for transcript_id in re.split(r"[,\s]+", parent.strip()):
                    if transcript_id:
                        exons_by_transcript[transcript_id].append((start_i, end_i))

    return genes, transcripts, exons_by_transcript, gene_to_transcripts


def gene_symbol_variants(symbol: str) -> set[str]:
    norm = normalize_gene_symbol(symbol)
    variants = {norm}
    stripped = re.sub(r"^[a-z]{2,4}(?:[-_])?", "", norm)
    variants.add(stripped)
    base = re.sub(r"^[A-Za-z]{2,4}[-_]", "", symbol or "")
    variants.add(normalize_gene_symbol(base))
    return {item for item in variants if item}


def is_excluded_gene(gene: dict) -> bool:
    text = f"{gene.get('gene_symbol', '')} {gene.get('gene_name', '')}"
    if normalize_gene_symbol(gene.get("gene_symbol", "")) in EXCLUDED_GENE_KEYWORDS:
        return True
    return bool(EXCLUDED_GENE_PATTERNS.search(text))


def build_gene_index(genes: dict[str, dict]) -> dict[str, list[str]]:
    index: dict[str, list[str]] = defaultdict(list)
    for gene_id, gene in genes.items():
        if is_excluded_gene(gene):
            continue

        keys = {normalize_gene_symbol(gene["gene_symbol"])}
        keys |= gene_symbol_variants(gene["gene_symbol"])
        keys |= product_search_keys(gene.get("gene_name", ""))

        for synonym in re.split(r"[,;]\s*", gene.get("gene_synonym", "")):
            synonym = synonym.strip()
            if synonym:
                keys.add(normalize_gene_symbol(synonym))
                keys |= gene_symbol_variants(synonym)

        for key in keys:
            if key:
                index[key].append(gene_id)

    return index


def product_search_keys(product: str) -> set[str]:
    keys: set[str] = set()
    if not product:
        return keys

    norm = normalize_gene_symbol(product)
    if norm:
        keys.add(norm)

    for token in re.split(r"[^a-zA-Z0-9]+", product.lower()):
        if not token:
            continue
        keys.add(normalize_gene_symbol(token))
        keys |= gene_symbol_variants(token)

    return {key for key in keys if key}


def build_transcript_product_index(
    transcripts: dict[str, dict],
    genes: dict[str, dict],
) -> dict[str, list[str]]:
    index: dict[str, list[str]] = defaultdict(list)
    for transcript in transcripts.values():
        parent = transcript.get("parent_gene")
        if not parent or parent not in genes or is_excluded_gene(genes[parent]):
            continue
        for key in product_search_keys(transcript.get("transcript_product", "")):
            index[key].append(parent)
    return index


def description_match_score(reference: dict[str, str], gene: dict) -> int:
    protein_name = (reference.get("protein_name") or "").lower()
    gene_text = f"{gene.get('gene_name', '')} {gene.get('gene_symbol', '')}".lower()
    if not protein_name or not gene_text:
        return 0

    score = 0
    if "gap junction" in protein_name and "gap junction" in gene_text:
        score += 2
    if "innexin" in protein_name and (
        "innexin" in gene_text or "inx" in gene_text or normalize_gene_symbol(gene["gene_symbol"]).startswith("inx")
    ):
        score += 2

    for token in ("alpha", "beta", "gamma", "delta"):
        if token in protein_name and token in gene_text:
            score += 1

    return score


def choose_best_gene_match(
    matches: list[str],
    genes: dict[str, dict],
    gene_to_transcripts: dict[str, list[str]],
) -> str:
    def rank(gene_id: str) -> tuple[int, int, int]:
        gene = genes[gene_id]
        biotype = (gene.get("gene_biotype") or "").lower()
        coding = 1 if biotype == "protein_coding" else 0
        has_tx = 1 if gene_to_transcripts.get(gene_id) else 0
        excluded = 0 if not is_excluded_gene(gene) else -5
        return (excluded, coding, has_tx)

    return sorted(matches, key=rank, reverse=True)[0]


def find_gene_for_reference(
    reference: dict[str, str],
    genes: dict[str, dict],
    gene_index: dict[str, list[str]],
    transcript_index: dict[str, list[str]],
    gene_to_transcripts: dict[str, list[str]],
) -> tuple[str | None, str]:
    for symbol in candidate_gene_symbols(reference):
        for variant in gene_symbol_variants(symbol):
            matches = [gene_id for gene_id in gene_index.get(variant, []) if gene_id in genes]
            if matches:
                return (
                    choose_best_gene_match(matches, genes, gene_to_transcripts),
                    f"gene_symbol:{symbol}",
                )

            product_matches = [
                gene_id for gene_id in transcript_index.get(variant, []) if gene_id in genes
            ]
            if product_matches:
                return (
                    choose_best_gene_match(product_matches, genes, gene_to_transcripts),
                    f"transcript_product:{symbol}",
                )

    protein_norm = normalize_gene_symbol(reference.get("protein_name", ""))
    if protein_norm:
        product_matches = [
            gene_id for gene_id in transcript_index.get(protein_norm, []) if gene_id in genes
        ]
        if product_matches:
            return (
                choose_best_gene_match(product_matches, genes, gene_to_transcripts),
                "transcript_product:protein_name",
            )

    scored: list[tuple[int, str]] = []
    for gene_id, gene in genes.items():
        if is_excluded_gene(gene):
            continue
        score = description_match_score(reference, gene)
        if score > 0:
            scored.append((score, gene_id))

    if scored:
        scored.sort(key=lambda item: (-item[0], item[1]))
        best_score = scored[0][0]
        best = [gene_id for score, gene_id in scored if score == best_score]
        if best:
            return (
                choose_best_gene_match(best, genes, gene_to_transcripts),
                "description_match",
            )

    return None, "no_gff_gene_match"


def compute_transcript_metrics(exons: list[tuple[int, int]]) -> dict[str, object]:
    exons_sorted = sorted(exons, key=lambda item: item[0])
    exon_lengths = [end - start + 1 for start, end in exons_sorted]
    intron_lengths = []

    for index in range(len(exons_sorted) - 1):
        intron_len = exons_sorted[index + 1][0] - exons_sorted[index][1] - 1
        if intron_len >= 0:
            intron_lengths.append(intron_len)

    transcript_span = (
        exons_sorted[-1][1] - exons_sorted[0][0] + 1 if exons_sorted else 0
    )

    return {
        "exon_count": len(exons_sorted),
        "mean_exon_length": mean(exon_lengths) if exon_lengths else 0,
        "total_exon_bases": sum(exon_lengths),
        "intron_count": len(intron_lengths),
        "mean_intron_length": mean(intron_lengths) if intron_lengths else 0,
        "transcript_span": transcript_span,
    }


def main(species_csv: str) -> None:
    species_rows = load_species_csv(Path(species_csv))

    detailed_rows: list[dict] = []
    transcript_rows: list[dict] = []
    gene_rows: list[dict] = []
    protein_rows: list[dict] = []
    not_found_rows: list[dict] = []

    for index, row in enumerate(species_rows, start=1):
        family = row["family"].lower()
        organism = row["organism"]
        print(f"\n[{index}/{len(species_rows)}] {organism} ({family})")

        references = load_reference_proteins(family, organism)
        if not references:
            not_found_rows.append(
                {
                    "family": family,
                    "organism": organism,
                    "fasta_path": "",
                    "uniprot_accession": "",
                    "gene_symbol": "",
                    "reason": "No reference FASTA files found",
                }
            )
            continue

        annotation_path = find_annotation_file(family, organism)
        if annotation_path is None:
            for reference in references:
                not_found_rows.append(
                    {
                        "family": family,
                        "organism": organism,
                        "fasta_path": reference["fasta_path"],
                        "uniprot_accession": reference["uniprot_accession"],
                        "gene_symbol": reference["gene_symbol"],
                        "reason": "No annotation file found",
                    }
                )
            continue

        genes, transcripts, exons_by_transcript, gene_to_transcripts = parse_annotation(
            annotation_path
        )
        gene_index = build_gene_index(genes)
        transcript_index = build_transcript_product_index(transcripts, genes)

        for reference in references:
            gene_id, match_method = find_gene_for_reference(
                reference,
                genes,
                gene_index,
                transcript_index,
                gene_to_transcripts,
            )

            protein_row = {
                "family": family,
                "organism": organism,
                "fasta_path": reference["fasta_path"],
                "uniprot_accession": reference["uniprot_accession"],
                "reference_gene_symbol": reference["gene_symbol"],
                "protein_name": reference["protein_name"],
                "matched_gene_id": gene_id or "",
                "match_method": match_method,
                "status": "matched" if gene_id else "unmatched",
            }
            protein_rows.append(protein_row)

            if not gene_id:
                not_found_rows.append(
                    {
                        "family": family,
                        "organism": organism,
                        "fasta_path": reference["fasta_path"],
                        "uniprot_accession": reference["uniprot_accession"],
                        "gene_symbol": reference["gene_symbol"],
                        "reason": match_method,
                    }
                )
                continue

            gene = genes[gene_id]
            transcript_ids = gene_to_transcripts.get(gene_id, [])
            exon_counts: list[int] = []

            if not transcript_ids:
                not_found_rows.append(
                    {
                        "family": family,
                        "organism": organism,
                        "fasta_path": reference["fasta_path"],
                        "uniprot_accession": reference["uniprot_accession"],
                        "gene_symbol": gene["gene_symbol"],
                        "reason": "Gene matched but no transcripts found",
                    }
                )

            for transcript_id in transcript_ids:
                transcript = transcripts[transcript_id]
                exons = sorted(
                    exons_by_transcript.get(transcript_id, []), key=lambda item: item[0]
                )
                metrics = compute_transcript_metrics(exons)
                exon_counts.append(metrics["exon_count"])

                for exon_number, (start, end) in enumerate(exons, start=1):
                    detailed_rows.append(
                        {
                            "family": family,
                            "organism": organism,
                            "fasta_path": reference["fasta_path"],
                            "uniprot_accession": reference["uniprot_accession"],
                            "reference_gene_symbol": reference["gene_symbol"],
                            "gene_id": gene_id,
                            "gene_symbol": gene["gene_symbol"],
                            "gene_name": gene["gene_name"],
                            "transcript_id": transcript_id,
                            "transcript_name": transcript["transcript_name"],
                            "seqid": gene["seqid"],
                            "strand": gene["strand"],
                            "exon_number": exon_number,
                            "exon_start": start,
                            "exon_end": end,
                            "exon_length": end - start + 1,
                        }
                    )

                transcript_rows.append(
                    {
                        "family": family,
                        "organism": organism,
                        "fasta_path": reference["fasta_path"],
                        "uniprot_accession": reference["uniprot_accession"],
                        "reference_gene_symbol": reference["gene_symbol"],
                        "gene_id": gene_id,
                        "gene_symbol": gene["gene_symbol"],
                        "gene_name": gene["gene_name"],
                        "transcript_id": transcript_id,
                        "transcript_name": transcript["transcript_name"],
                        "seqid": gene["seqid"],
                        "strand": gene["strand"],
                        "transcript_start": transcript["start"],
                        "transcript_end": transcript["end"],
                        "transcript_span": metrics["transcript_span"],
                        "exon_count": metrics["exon_count"],
                        "mean_exon_length": metrics["mean_exon_length"],
                        "total_exon_bases": metrics["total_exon_bases"],
                        "intron_count": metrics["intron_count"],
                        "mean_intron_length": metrics["mean_intron_length"],
                    }
                )

            gene_rows.append(
                {
                    "family": family,
                    "organism": organism,
                    "fasta_path": reference["fasta_path"],
                    "uniprot_accession": reference["uniprot_accession"],
                    "reference_gene_symbol": reference["gene_symbol"],
                    "gene_id": gene_id,
                    "gene_symbol": gene["gene_symbol"],
                    "gene_name": gene["gene_name"],
                    "seqid": gene["seqid"],
                    "strand": gene["strand"],
                    "gene_start": gene["start"],
                    "gene_end": gene["end"],
                    "gene_span": gene["end"] - gene["start"] + 1,
                    "transcript_count": len(transcript_ids),
                    "mean_exon_count": mean(exon_counts) if exon_counts else 0,
                    "max_exon_count": max(exon_counts) if exon_counts else 0,
                    "min_exon_count": min(exon_counts) if exon_counts else 0,
                }
            )

    outputs = {
        "exon_coordinates.csv": detailed_rows,
        "transcript_summary.csv": transcript_rows,
        "gene_summary.csv": gene_rows,
        "protein_mapping.csv": protein_rows,
        "not_found_or_unmatched.csv": not_found_rows,
    }

    for filename, rows in outputs.items():
        pd.DataFrame(rows).to_csv(RESULTS_DIR / filename, index=False)
        print(f"Saved {filename} ({len(rows)} rows)")

    matched = sum(1 for row in protein_rows if row["status"] == "matched")
    print(f"\nMatched {matched}/{len(protein_rows)} reference proteins to GFF genes")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python extract_exon_structures.py <species.csv>")
        sys.exit(1)
    main(sys.argv[1])
