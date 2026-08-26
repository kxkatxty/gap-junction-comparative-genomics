from __future__ import annotations

import argparse
import concurrent.futures
import csv
import gzip
import json
import os
import re
import shutil
import subprocess
import sys
import tempfile
from collections import defaultdict
from dataclasses import dataclass, field
from pathlib import Path
from urllib.parse import urlparse

from download_gff3_from_species_list import (
    download_file,
    fetch_assembly_candidates,
    rank_assemblies,
    sanitize_filename,
)


PROJECT_ROOT = Path("project")

FAMILY_DEFAULTS: dict[str, dict[str, Path | str]] = {
    "innexin": {
        "references_dir": PROJECT_ROOT / "data" / "references" / "innexins",
        "annotation_dir": PROJECT_ROOT / "data" / "annotations" / "innexin",
        "genome_dir": PROJECT_ROOT / "data" / "genomes" / "innexin_discovery",
        "results_dir": PROJECT_ROOT / "results" / "innexin_discovery",
        "species_txt": PROJECT_ROOT / "metadata" / "not_annotated_innexin.txt",
        "trusted_fasta_name": "trusted_innexins.fasta",
        "candidate_fasta_name": "candidate_innexins.fasta",
        "candidate_gff3_name": "candidate_innexins.gff3",
        "plot_name": "candidate_innexin_exon_map.png",
    },
    "connexin": {
        "references_dir": PROJECT_ROOT / "data" / "references" / "connexins",
        "annotation_dir": PROJECT_ROOT / "data" / "annotations" / "connexin",
        "genome_dir": PROJECT_ROOT / "data" / "genomes" / "innexin_discovery",
        "results_dir": PROJECT_ROOT / "results" / "connexin_discovery",
        "species_txt": PROJECT_ROOT / "metadata" / "not_annotated_connexin.txt",
        "trusted_fasta_name": "trusted_connexins.fasta",
        "candidate_fasta_name": "candidate_connexins.fasta",
        "candidate_gff3_name": "candidate_connexins.gff3",
        "plot_name": "candidate_connexin_exon_map.png",
    },
}

ACTIVE_FAMILY = "innexin"
REFERENCES_DIR = FAMILY_DEFAULTS["innexin"]["references_dir"]  # type: ignore[assignment]
ANNOTATION_DIR = FAMILY_DEFAULTS["innexin"]["annotation_dir"]  # type: ignore[assignment]
GENOME_DIR = FAMILY_DEFAULTS["innexin"]["genome_dir"]  # type: ignore[assignment]
RESULTS_DIR = FAMILY_DEFAULTS["innexin"]["results_dir"]  # type: ignore[assignment]
DEFAULT_SPECIES_TXT = FAMILY_DEFAULTS["innexin"]["species_txt"]  # type: ignore[assignment]
TRUSTED_FASTA_NAME = "trusted_innexins.fasta"
CANDIDATE_FASTA_NAME = "candidate_innexins.fasta"
CANDIDATE_GFF3_NAME = "candidate_innexins.gff3"
PLOT_NAME = "candidate_innexin_exon_map.png"


def configure_family(family: str) -> None:
    global ACTIVE_FAMILY, REFERENCES_DIR, ANNOTATION_DIR, GENOME_DIR, RESULTS_DIR
    global DEFAULT_SPECIES_TXT, TRUSTED_FASTA_NAME, CANDIDATE_FASTA_NAME
    global CANDIDATE_GFF3_NAME, PLOT_NAME, RANK_HIGH, RANK_POSSIBLE, RANK_WEAK, RANK_REJECTED

    key = family.lower().strip()
    if key not in FAMILY_DEFAULTS:
        raise ValueError(f"Unknown family {family!r}; choose from: {', '.join(FAMILY_DEFAULTS)}")

    profile = FAMILY_DEFAULTS[key]
    ACTIVE_FAMILY = key
    REFERENCES_DIR = profile["references_dir"]  # type: ignore[assignment]
    ANNOTATION_DIR = profile["annotation_dir"]  # type: ignore[assignment]
    GENOME_DIR = profile["genome_dir"]  # type: ignore[assignment]
    RESULTS_DIR = profile["results_dir"]  # type: ignore[assignment]
    DEFAULT_SPECIES_TXT = profile["species_txt"]  # type: ignore[assignment]
    TRUSTED_FASTA_NAME = str(profile["trusted_fasta_name"])
    CANDIDATE_FASTA_NAME = str(profile["candidate_fasta_name"])
    CANDIDATE_GFF3_NAME = str(profile["candidate_gff3_name"])
    PLOT_NAME = str(profile["plot_name"])
    RANK_HIGH = f"high_confidence_{key}_candidate"
    RANK_POSSIBLE = f"possible_{key}_candidate"
    RANK_WEAK = f"weak_manual_review"
    RANK_REJECTED = "rejected_false_positive"

CLUSTER_MAX_GAP_BP = 8_000
FLANK_BP = 2_000

# miniprot genome search (protein queries vs nucleotide genome; splice-aware)
MINIPROT_SPLICE_MODEL = 2  # 2 = insect/vertebrate; good for arthropods/annelids
MINIPROT_MIN_SCORE = 30.0  # miniprot alignment score floor
CHUNKED_GENOME_BYTES = 4 * 1024**3  # per-contig miniprot above this genome size
CHUNKED_MIN_CONTIG_BP = 1_000_000  # large scaffolds/chromosomes only in chunked mode
CHUNKED_MAX_WORKERS = 3  # parallel miniprot jobs per large genome

# Post-search HSP filters (relaxed for distant/partial innexin homologs)
MIN_HSP_LENGTH_BP = 60  # ≥ ~20 codons on the genome
MIN_HSP_QUERY_AA = 20  # ~partial TM-region alignment
MIN_HSP_PIDENT = 18.0  # allow more divergent miniprot hits
MIN_HSP_BITSCORE = 25.0
MIN_HSP_SHORT_QUERY_AA = 35
MIN_HSP_SHORT_PIDENT = 22.0

# Candidate ORF / validation
MIN_CANDIDATE_ORF_AA = 70  # allow shorter partial innexins
BLASTP_EVALUE = "1e-5"

# MMseqs2 protein-vs-protein validation only (candidate ranking)
MMSEQS_SEARCH_TYPE_BLASTP = 1

SEARCH_THREADS = max(1, os.cpu_count() or 1)
MMSEQS_BLASTP_FORMAT = "query,target,pident,qcov,bits"

MINIPROT_TARGET_RE = re.compile(r"Target=([^\s;]+)\s+(\d+)\s+(\d+)")
MINIPROT_IDENTITY_RE = re.compile(r"Identity=([\d.]+)")

INNEXIN_MIN_AA = 150
INNEXIN_MAX_AA = 750
MIN_TM_HELICES = 2  # innexins have 4 TMs; accept 2+ for partial candidates
MIN_REFERENCE_IDENTITY = 15.0
MIN_INNEXIN_REJECT_AA = 120  # hard reject below this length in scoring

RANK_HIGH = "high_confidence_innexin_candidate"
RANK_POSSIBLE = "possible_innexin_candidate"
RANK_WEAK = "weak_manual_review"
RANK_REJECTED = "rejected_false_positive"

FALSE_POSITIVE_KEYWORDS: tuple[str, ...] = (
    "kinase",
    "receptor tyrosine",
    "serine/threonine-protein kinase",
    "tyrosine-protein kinase",
    "toll-like",
    "disease resistance",
    "nbs-lrr",
    "leucine-rich repeat",
    "immunoglobulin",
    "histone",
    "helicase",
    "wd40",
    "ankyrin repeat",
    "g-protein coupled",
    "gpcr",
    "sodium channel",
    "potassium channel",
    "calcium channel",
    "abc transporter",
    "mhc class",
)

KYTE_DOOLITTLE = {
    "A": 1.8,
    "C": 2.5,
    "D": -3.5,
    "E": -3.5,
    "F": 2.8,
    "G": -0.4,
    "H": -3.2,
    "I": 4.5,
    "K": -3.9,
    "L": 3.8,
    "M": 1.9,
    "N": -3.5,
    "P": -1.6,
    "Q": -3.5,
    "R": -4.5,
    "S": -0.8,
    "T": -0.7,
    "V": 4.2,
    "W": -0.9,
    "Y": -1.3,
}

GENETIC_CODE = {
    "TTT": "F",
    "TTC": "F",
    "TTA": "L",
    "TTG": "L",
    "TCT": "S",
    "TCC": "S",
    "TCA": "S",
    "TCG": "S",
    "TAT": "Y",
    "TAC": "Y",
    "TAA": "*",
    "TAG": "*",
    "TGT": "C",
    "TGC": "C",
    "TGA": "*",
    "TGG": "W",
    "CTT": "L",
    "CTC": "L",
    "CTA": "L",
    "CTG": "L",
    "CCT": "P",
    "CCC": "P",
    "CCA": "P",
    "CCG": "P",
    "CAT": "H",
    "CAC": "H",
    "CAA": "Q",
    "CAG": "Q",
    "CGT": "R",
    "CGC": "R",
    "CGA": "R",
    "CGG": "R",
    "ATT": "I",
    "ATC": "I",
    "ATA": "I",
    "ATG": "M",
    "ACT": "T",
    "ACC": "T",
    "ACA": "T",
    "ACG": "T",
    "AAT": "N",
    "AAC": "N",
    "AAA": "K",
    "AAG": "K",
    "AGT": "S",
    "AGC": "S",
    "AGA": "R",
    "AGG": "R",
    "GTT": "V",
    "GTC": "V",
    "GTA": "V",
    "GTG": "V",
    "GCT": "A",
    "GCC": "A",
    "GCA": "A",
    "GCG": "A",
    "GAT": "D",
    "GAC": "D",
    "GAA": "E",
    "GAG": "E",
    "GGT": "G",
    "GGC": "G",
    "GGA": "G",
    "GGG": "G",
}


@dataclass
class BlastHsp:
    query_id: str
    target_id: str
    query_start: int
    query_end: int
    target_start: int
    target_end: int
    evalue: float
    bitscore: float
    identity: float
    frame: int
    qseq: str = ""
    tseq: str = ""


@dataclass
class ExonModel:
    exon_number: int
    start: int
    end: int
    phase: int = 0


@dataclass
class CandidateLocus:
    organism: str
    candidate_id: str
    seqid: str
    strand: str
    locus_start: int
    locus_end: int
    window_start: int
    window_end: int
    exons: list[ExonModel] = field(default_factory=list)
    protein_sequence: str = ""
    cds_sequence: str = ""
    protein_length: int = 0
    exon_count: int = 0
    tm_helix_count: int = 0
    reference_identity: float = 0.0
    reference_coverage: float = 0.0
    best_reference_hit: str = ""
    reference_hit_count: int = 0
    completeness: str = ""
    validation_flags: list[str] = field(default_factory=list)
    rank_category: str = ""
    rank_score: float = 0.0
    rejection_reason: str = ""
    hsps: list[BlastHsp] = field(default_factory=list)


def normalize_species_name(name: str) -> str:
    return sanitize_filename(name.strip())


def read_species_txt(path: Path) -> list[str]:
    species: list[str] = []
    for line in path.read_text(encoding="utf-8").splitlines():
        line = line.strip()
        if line and not line.startswith("#"):
            species.append(line)
    return species


def is_gzip_file(path: Path) -> bool:
    with path.open("rb") as handle:
        return handle.read(2) == b"\x1f\x8b"


def genome_fasta_is_valid(path: Path) -> bool:
    if not path.exists() or path.stat().st_size < 1000:
        return False
    try:
        if is_gzip_file(path):
            with gzip.open(path, "rb") as handle:
                while handle.read(1024 * 1024):
                    pass
        else:
            with path.open("rb") as handle:
                handle.read(1024)
        return True
    except Exception:
        return False


def decompress_if_gzip(path: Path) -> Path:
    if not is_gzip_file(path):
        return path
    out_path = path.with_suffix("") if path.suffix == ".gz" else path.with_name(path.name + ".decompressed.fna")
    if out_path.exists() and out_path.stat().st_size > 0 and genome_fasta_is_valid(out_path):
        return out_path
    if out_path.exists():
        out_path.unlink(missing_ok=True)
    with gzip.open(path, "rb") as src, out_path.open("wb") as dst:
        shutil.copyfileobj(src, dst)
    return out_path


def clear_genome_files(species_dir: Path) -> None:
    if not species_dir.exists():
        return
    for pattern in ("*_genomic.fna", "*_genomic.fna.gz", "*_genomic.fna.decompressed.fna", "*.fai"):
        for path in species_dir.glob(pattern):
            path.unlink(missing_ok=True)


def ensure_faidx(genome_path: Path) -> None:
    index_path = Path(f"{genome_path}.fai")
    if index_path.exists():
        return
    subprocess.run(["samtools", "faidx", str(genome_path)], check=True, capture_output=True)


def contig_length(genome_path: Path, seqid: str) -> int:
    ensure_faidx(genome_path)
    with Path(f"{genome_path}.fai").open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if parts and parts[0] == seqid:
                return int(parts[1])
    return 0


def extract_genome_region(genome_path: Path, seqid: str, start: int, end: int) -> str:
    ensure_faidx(genome_path)
    result = subprocess.run(
        ["samtools", "faidx", str(genome_path), f"{seqid}:{start}-{end}"],
        check=True,
        capture_output=True,
        text=True,
    )
    return "".join(
        line.strip()
        for line in result.stdout.splitlines()
        if line and not line.startswith(">")
    )


def read_fasta_records(path: Path) -> dict[str, str]:
    records: dict[str, str] = {}
    current_id = ""
    chunks: list[str] = []
    use_gzip = str(path).endswith(".gz") or is_gzip_file(path)
    opener = gzip.open if use_gzip else open
    with opener(path, "rt", encoding="utf-8") as handle:  # type: ignore[arg-type]
        for line in handle:
            if line.startswith(">"):
                if current_id:
                    records[current_id] = "".join(chunks)
                current_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.strip())
        if current_id:
            records[current_id] = "".join(chunks)
    return records


def write_fasta_records(records: dict[str, str], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for seq_id, sequence in records.items():
            handle.write(f">{seq_id}\n")
            for index in range(0, len(sequence), 80):
                handle.write(sequence[index : index + 80] + "\n")


def collect_reference_sequences(reference_dir: Path) -> tuple[Path, list[dict[str, str]]]:
    fasta_files = sorted(reference_dir.rglob("*.fasta"))
    if not fasta_files:
        raise FileNotFoundError(
            f"No {ACTIVE_FAMILY} FASTA files found in {reference_dir}"
        )

    combined_records: dict[str, str] = {}
    metadata: list[dict[str, str]] = []
    for fasta_path in fasta_files:
        for header_id, sequence in read_fasta_records(fasta_path).items():
            match = re.match(r"^(.+?)__([A-Z0-9]+)\.fasta$", fasta_path.name, re.I)
            gene_symbol = match.group(1) if match else fasta_path.stem
            accession = match.group(2).upper() if match else ""
            record_id = f"{gene_symbol}__{accession}" if accession else gene_symbol
            combined_records[record_id] = sequence
            metadata.append(
                {
                    "record_id": record_id,
                    "gene_symbol": gene_symbol,
                    "accession": accession,
                    "source_species": fasta_path.parent.name,
                    "fasta_path": str(fasta_path),
                }
            )

    out_dir = RESULTS_DIR / "_cache"
    out_dir.mkdir(parents=True, exist_ok=True)
    combined_path = out_dir / TRUSTED_FASTA_NAME
    write_fasta_records(combined_records, combined_path)
    return combined_path, metadata


def collect_reference_innexins(reference_dir: Path) -> tuple[Path, list[dict[str, str]]]:
    return collect_reference_sequences(reference_dir)


def tool_exists(name: str) -> bool:
    return shutil.which(name) is not None


def require_discovery_tools() -> None:
    missing = [tool for tool in ("miniprot", "mmseqs", "samtools") if not tool_exists(tool)]
    if missing:
        raise RuntimeError(
            "Discovery tools missing: "
            + ", ".join(missing)
            + ". Install with: conda install -c bioconda miniprot mmseqs2 samtools"
        )


def genome_fna_url_from_assembly(assembly: dict[str, str]) -> str | None:
    for source_key in ("FtpPath_RefSeq", "FtpPath_GenBank"):
        ftp_path = assembly.get(source_key, "")
        if not ftp_path:
            continue
        assembly_base = ftp_path.rstrip("/").split("/")[-1]
        https_path = ftp_path.replace("ftp://", "https://")
        return f"{https_path}/{assembly_base}_genomic.fna.gz"
    return None


def ensure_genome_fasta(organism: str, genome_dir: Path, *, force_redownload: bool = False) -> Path:
    species_dir = genome_dir / normalize_species_name(organism)
    species_dir.mkdir(parents=True, exist_ok=True)

    if force_redownload:
        clear_genome_files(species_dir)

    existing = sorted(species_dir.glob("*_genomic.fna")) + sorted(
        species_dir.glob("*_genomic.fna.gz")
    )
    for genome_path in existing:
        if not genome_fasta_is_valid(genome_path):
            clear_genome_files(species_dir)
            break
        try:
            return decompress_if_gzip(genome_path)
        except Exception:
            clear_genome_files(species_dir)
            break

    candidates = rank_assemblies(fetch_assembly_candidates(organism))
    if not candidates:
        raise FileNotFoundError(f"No genome assembly found for {organism}")

    errors: list[str] = []
    for assembly in candidates[:8]:
        url = genome_fna_url_from_assembly(assembly)
        if not url:
            continue
        assembly_base = Path(urlparse(url).path).name.replace(".gz", "")
        out_path = species_dir / assembly_base
        out_path.unlink(missing_ok=True)
        try:
            download_file(url, out_path)
            return decompress_if_gzip(out_path)
        except Exception as exc:  # noqa: BLE001
            out_path.unlink(missing_ok=True)
            errors.append(f"{assembly.get('AssemblyAccession')}: {exc}")

    raise FileNotFoundError(
        f"Could not download genome for {organism}. Errors: {'; '.join(errors)}"
    )


def run_mmseqs_easy_search(
    query_fasta: Path,
    target_fasta: Path,
    out_m8: Path,
    tmp_dir: Path,
    *,
    search_type: int,
    evalue: str,
    max_seqs: int,
    format_output: str,
    sensitivity: float = 5.7,
    num_threads: int | None = None,
) -> None:
    threads = num_threads if num_threads is not None else SEARCH_THREADS
    tmp_dir.mkdir(parents=True, exist_ok=True)
    out_m8.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mmseqs",
        "easy-search",
        str(query_fasta),
        str(target_fasta),
        str(out_m8),
        str(tmp_dir),
        "--search-type",
        str(search_type),
        "-e",
        evalue,
        "--max-seqs",
        str(max_seqs),
        "-s",
        str(sensitivity),
        "--threads",
        str(threads),
        "--format-output",
        format_output,
    ]
    subprocess.run(cmd, check=True, capture_output=True, text=True)


def list_genome_contigs(genome_fasta: Path) -> list[tuple[str, int]]:
    ensure_faidx(genome_fasta)
    contigs: list[tuple[str, int]] = []
    with Path(f"{genome_fasta}.fai").open("r", encoding="utf-8") as handle:
        for line in handle:
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 2:
                contigs.append((parts[0], int(parts[1])))
    contigs.sort(key=lambda item: item[1], reverse=True)
    return contigs


def run_miniprot_genome_search(
    genome_fasta: Path,
    query_fasta: Path,
    out_gff: Path,
    *,
    num_threads: int | None = None,
    splice_model: int = MINIPROT_SPLICE_MODEL,
) -> None:
    threads = num_threads if num_threads is not None else SEARCH_THREADS
    out_gff.parent.mkdir(parents=True, exist_ok=True)
    if genome_fasta.stat().st_size >= CHUNKED_GENOME_BYTES:
        run_miniprot_genome_search_chunked(
            genome_fasta,
            query_fasta,
            out_gff,
            num_threads=threads,
            splice_model=splice_model,
        )
        return
    cmd = [
        "miniprot",
        f"-t{threads}",
        f"-j{splice_model}",
        "--gff",
        str(genome_fasta),
        str(query_fasta),
    ]
    with out_gff.open("w", encoding="utf-8") as handle:
        subprocess.run(cmd, check=True, stdout=handle, stderr=subprocess.PIPE, text=True)


def _run_single_miniprot_chunk(
    genome_fasta: Path,
    seqid: str,
    chunk_dir: Path,
    query_fasta: Path,
    threads: int,
    splice_model: int,
) -> str:
    chunk_fasta = chunk_dir / f"{seqid}.fa"
    chunk_gff = chunk_dir / f"{seqid}.gff"
    if chunk_gff.exists() and chunk_gff.stat().st_size > 0:
        return chunk_gff.read_text(encoding="utf-8")

    if not chunk_fasta.exists() or chunk_fasta.stat().st_size == 0:
        result = subprocess.run(
            ["samtools", "faidx", str(genome_fasta), seqid],
            check=True,
            capture_output=True,
            text=True,
        )
        chunk_fasta.write_text(result.stdout, encoding="utf-8")

    cmd = [
        "miniprot",
        f"-t{threads}",
        f"-j{splice_model}",
        "--gff",
        str(chunk_fasta),
        str(query_fasta),
    ]
    with chunk_gff.open("w", encoding="utf-8") as handle:
        subprocess.run(cmd, check=True, stdout=handle, stderr=subprocess.PIPE, text=True)
    if chunk_gff.stat().st_size > 0:
        return chunk_gff.read_text(encoding="utf-8")
    return ""


def run_miniprot_genome_search_chunked(
    genome_fasta: Path,
    query_fasta: Path,
    out_gff: Path,
    *,
    num_threads: int | None = None,
    splice_model: int = MINIPROT_SPLICE_MODEL,
    min_contig_bp: int = CHUNKED_MIN_CONTIG_BP,
    max_workers: int = CHUNKED_MAX_WORKERS,
) -> None:
    per_chunk_threads = max(1, min(num_threads or SEARCH_THREADS, 2))
    contigs = [
        (seqid, length)
        for seqid, length in list_genome_contigs(genome_fasta)
        if length >= min_contig_bp
    ]
    if not contigs:
        raise RuntimeError(f"No contigs >= {min_contig_bp} bp in {genome_fasta}")

    chunk_dir = out_gff.parent / "miniprot_chunks"
    chunk_dir.mkdir(parents=True, exist_ok=True)
    workers = max(1, min(max_workers, len(contigs)))
    parts: list[str] = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        futures = [
            executor.submit(
                _run_single_miniprot_chunk,
                genome_fasta,
                seqid,
                chunk_dir,
                query_fasta,
                per_chunk_threads,
                splice_model,
            )
            for seqid, _length in contigs
        ]
        for future in concurrent.futures.as_completed(futures):
            text = future.result()
            if text:
                parts.append(text)
    out_gff.write_text("".join(parts), encoding="utf-8")


def parse_miniprot_gff(path: Path) -> list[BlastHsp]:
    hsps: list[BlastHsp] = []
    if not path.exists() or path.stat().st_size == 0:
        return hsps

    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9 or parts[2] != "mRNA":
                continue

            target_match = MINIPROT_TARGET_RE.search(parts[8])
            if not target_match:
                continue

            identity_match = MINIPROT_IDENTITY_RE.search(parts[8])
            identity = float(identity_match.group(1)) * 100.0 if identity_match else 0.0
            score = float(parts[5]) if parts[5] not in {"", "."} else 0.0
            strand = parts[6]
            target_start = int(parts[3])
            target_end = int(parts[4])
            frame = -1 if strand == "-" else 1

            hsps.append(
                BlastHsp(
                    query_id=target_match.group(1),
                    target_id=parts[0],
                    query_start=int(target_match.group(2)),
                    query_end=int(target_match.group(3)),
                    target_start=min(target_start, target_end),
                    target_end=max(target_start, target_end),
                    evalue=1e-20,
                    bitscore=score,
                    identity=identity,
                    frame=frame,
                )
            )
    return hsps


def hsp_target_span_bp(hsp: BlastHsp) -> int:
    return hsp.target_end - hsp.target_start + 1


def hsp_query_span_aa(hsp: BlastHsp) -> int:
    return abs(hsp.query_end - hsp.query_start) + 1


def passes_hsp_quality_filters(hsp: BlastHsp) -> bool:
    query_aa = hsp_query_span_aa(hsp)
    if hsp_target_span_bp(hsp) < MIN_HSP_LENGTH_BP:
        return False
    if query_aa < MIN_HSP_QUERY_AA:
        return False
    if hsp.bitscore < max(MIN_HSP_BITSCORE, MINIPROT_MIN_SCORE):
        return False
    if hsp.identity < MIN_HSP_PIDENT:
        return False
    if query_aa < MIN_HSP_SHORT_QUERY_AA and hsp.identity < MIN_HSP_SHORT_PIDENT:
        return False
    return True


def filter_hsps(hsps: list[BlastHsp]) -> list[BlastHsp]:
    return [hsp for hsp in hsps if passes_hsp_quality_filters(hsp)]


def hsp_strand(hsp: BlastHsp) -> str:
    return "-" if hsp.frame < 0 else "+"


def cluster_hsps_into_loci(
    hsps: list[BlastHsp],
    *,
    max_gap_bp: int = CLUSTER_MAX_GAP_BP,
) -> list[list[BlastHsp]]:
    if not hsps:
        return []

    grouped: dict[tuple[str, str], list[BlastHsp]] = defaultdict(list)
    for hsp in hsps:
        if not passes_hsp_quality_filters(hsp):
            continue
        grouped[(hsp.target_id, hsp_strand(hsp))].append(hsp)

    loci: list[list[BlastHsp]] = []
    for (_seqid, _strand), items in grouped.items():
        items.sort(key=lambda item: item.target_start)
        current = [items[0]]
        for hsp in items[1:]:
            prev = current[-1]
            if hsp.target_start <= prev.target_end + max_gap_bp:
                current.append(hsp)
            else:
                loci.append(current)
                current = [hsp]
        loci.append(current)
    return loci


def find_annotation_gff(organism: str) -> Path | None:
    species_dir = ANNOTATION_DIR / normalize_species_name(organism)
    if not species_dir.exists():
        return None
    candidates = sorted(species_dir.glob("*.gff")) + sorted(species_dir.glob("*.gff3"))
    candidates += sorted(species_dir.glob("*.gff.gz")) + sorted(species_dir.glob("*.gff3.gz"))
    return candidates[0] if candidates else None


def parse_gff_attributes(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if "=" in attr_str:
        for field in attr_str.split(";"):
            if "=" in field:
                key, value = field.split("=", 1)
                attrs[key.strip()] = value.strip().strip('"')
    else:
        for match in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
            attrs[match.group(1)] = match.group(2)
    return attrs


def gff_features_overlapping(
    gff_path: Path, seqid: str, start: int, end: int
) -> list[dict[str, str]]:
    overlaps: list[dict[str, str]] = []
    opener = gzip.open if str(gff_path).endswith(".gz") else open
    with opener(gff_path, "rt", encoding="utf-8") as handle:  # type: ignore[arg-type]
        for line in handle:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            chrom, _source, feature, fstart, fend, _score, _strand, _phase, attrs = parts
            if chrom != seqid:
                continue
            fs = int(fstart)
            fe = int(fend)
            if fe < start or fs > end:
                continue
            overlaps.append(
                {
                    "feature": feature,
                    "start": str(fs),
                    "end": str(fe),
                    "attributes": attrs,
                }
            )
    return overlaps


def translate_dna(sequence: str, frame: int = 0) -> str:
    seq = sequence.upper()
    if frame < 0:
        comp = str.maketrans("ACGT", "TGCA")
        seq = seq.translate(comp)[::-1]
        frame = abs(frame) - 1
    else:
        frame = max(frame - 1, 0)

    protein: list[str] = []
    for index in range(frame, len(seq) - 2, 3):
        codon = seq[index : index + 3]
        if len(codon) < 3 or "N" in codon:
            continue
        aa = GENETIC_CODE.get(codon, "X")
        if aa == "*":
            break
        protein.append(aa)
    return "".join(protein)


def build_exon_model_from_hsps(hsps: list[BlastHsp]) -> list[ExonModel]:
    exon_intervals: list[tuple[int, int]] = []
    for hsp in sorted(hsps, key=lambda item: item.target_start):
        if not exon_intervals:
            exon_intervals.append((hsp.target_start, hsp.target_end))
            continue
        last_start, last_end = exon_intervals[-1]
        if hsp.target_start <= last_end + 50:
            exon_intervals[-1] = (last_start, max(last_end, hsp.target_end))
        else:
            exon_intervals.append((hsp.target_start, hsp.target_end))

    exons: list[ExonModel] = []
    for index, (start, end) in enumerate(exon_intervals, start=1):
        exons.append(ExonModel(exon_number=index, start=start, end=end))
    return exons


def longest_orf_protein(window_seq: str, hsps: list[BlastHsp]) -> tuple[str, str, str]:
    best_protein = ""
    best_cds = ""
    best_frame = 0
    for frame in (1, 2, 3, -1, -2, -3):
        protein = translate_dna(window_seq, frame)
        if len(protein) <= len(best_protein):
            continue
        overlap = False
        for hsp in hsps:
            # keep frame that yields longest protein overlapping HSP span
            overlap = True
        if overlap:
            best_protein = protein
            best_cds = window_seq
            best_frame = frame
    completeness = "complete" if best_protein.startswith("M") and "*" not in best_protein else "partial"
    return best_protein, best_cds, completeness


def predict_tm_helices(protein: str, window: int = 19, threshold: float = 1.15) -> int:
    if len(protein) < window:
        return 0
    scores = [
        sum(KYTE_DOOLITTLE.get(aa, 0.0) for aa in protein[i : i + window]) / window
        for i in range(len(protein) - window + 1)
    ]
    helices = 0
    in_helix = False
    for score in scores:
        if score >= threshold and not in_helix:
            helices += 1
            in_helix = True
        elif score < threshold:
            in_helix = False
    return helices


def mmseqs_best_hit(
    query_fasta: Path,
    subject_fasta: Path,
    work_dir: Path,
    *,
    num_threads: int | None = None,
) -> tuple[float, float, str]:
    out_tsv = work_dir / f"mmseqs_{query_fasta.stem}.m8"
    tmp_dir = work_dir / f"mmseqs_{query_fasta.stem}_tmp"
    try:
        run_mmseqs_easy_search(
            query_fasta,
            subject_fasta,
            out_tsv,
            tmp_dir,
            search_type=MMSEQS_SEARCH_TYPE_BLASTP,
            evalue=BLASTP_EVALUE,
            max_seqs=1,
            format_output=MMSEQS_BLASTP_FORMAT,
            num_threads=num_threads,
        )
    except subprocess.CalledProcessError:
        return 0.0, 0.0, ""
    if not out_tsv.exists() or out_tsv.stat().st_size == 0:
        return 0.0, 0.0, ""
    parts = out_tsv.read_text(encoding="utf-8").splitlines()[0].split("\t")
    return float(parts[2]), float(parts[3]), parts[1]


def contains_false_positive_keyword(text: str) -> str | None:
    lowered = text.lower()
    for keyword in FALSE_POSITIVE_KEYWORDS:
        if keyword in lowered:
            return keyword
    return None


def score_candidate(candidate: CandidateLocus) -> tuple[float, str, str]:
    score = 0.0
    flags: list[str] = []
    rejection = ""

    if INNEXIN_MIN_AA <= candidate.protein_length <= INNEXIN_MAX_AA:
        score += 20
        flags.append("plausible_length")
    elif candidate.protein_length < INNEXIN_MIN_AA:
        flags.append("too_short")
        rejection = "protein too short"
    else:
        flags.append("long_protein")

    if candidate.tm_helix_count >= MIN_TM_HELICES:
        score += 25
        flags.append("tm_helices_ok")
    else:
        flags.append("few_tm_helices")

    if candidate.reference_identity >= MIN_REFERENCE_IDENTITY:
        score += 25
        flags.append("reference_similarity_ok")
    else:
        flags.append("low_reference_similarity")

    if candidate.exon_count >= 1:
        score += 10
        flags.append("gene_like_exons")

    if candidate.reference_hit_count >= 2:
        score += 10
        flags.append("multiple_reference_hits")

    if candidate.completeness == "complete":
        score += 10
        flags.append("complete_orf")

    fp_text = f"{candidate.best_reference_hit} {candidate.protein_sequence[:80]}"
    fp_keyword = contains_false_positive_keyword(fp_text)
    if not fp_keyword and ACTIVE_FAMILY == "connexin" and re.search(
        r"\binnexin\b", fp_text.lower()
    ):
        fp_keyword = "innexin"
    if fp_keyword:
        score -= 40
        flags.append(f"false_positive_keyword:{fp_keyword}")
        rejection = rejection or f"false positive keyword: {fp_keyword}"

    if candidate.protein_length < MIN_INNEXIN_REJECT_AA:
        score -= 30
        rejection = rejection or f"too short for {ACTIVE_FAMILY}"

    if candidate.tm_helix_count == 0:
        score -= 15
        if candidate.protein_length < INNEXIN_MIN_AA:
            rejection = rejection or "no transmembrane helices detected"

    if score >= 60 and not rejection:
        category = RANK_HIGH
    elif score >= 35 and not rejection:
        category = RANK_POSSIBLE
    elif score >= 15:
        category = RANK_WEAK
    else:
        category = RANK_REJECTED
        rejection = rejection or "low composite score"

    candidate.validation_flags = flags
    candidate.rank_score = round(score, 2)
    candidate.rank_category = category
    candidate.rejection_reason = rejection if category == RANK_REJECTED else ""
    return score, category, rejection


def build_candidate_locus(
    organism: str,
    locus_index: int,
    hsps: list[BlastHsp],
    genome_path: Path,
    *,
    flank_bp: int = FLANK_BP,
    gff_path: Path | None = None,
    work_dir: Path | None = None,
) -> CandidateLocus | None:
    if not hsps:
        return None

    seqid = hsps[0].target_id
    strand = hsp_strand(hsps[0])
    locus_start = min(hsp.target_start for hsp in hsps)
    locus_end = max(hsp.target_end for hsp in hsps)
    window_start = max(1, locus_start - flank_bp)
    contig_len = contig_length(genome_path, seqid)
    if contig_len == 0:
        return None
    window_end = min(contig_len, locus_end + flank_bp)
    window_seq = extract_genome_region(genome_path, seqid, window_start, window_end)
    if strand == "-":
        comp = str.maketrans("ACGTacgt", "TGCAtgca")
        window_seq = window_seq.translate(comp)[::-1]

    exons = build_exon_model_from_hsps(hsps)
    protein, cds, completeness = longest_orf_protein(window_seq, hsps)
    if len(protein) < MIN_CANDIDATE_ORF_AA:
        return None

    candidate = CandidateLocus(
        organism=organism,
        candidate_id=f"candidate_inx_{locus_index:03d}",
        seqid=seqid,
        strand=strand,
        locus_start=locus_start,
        locus_end=locus_end,
        window_start=window_start,
        window_end=window_end,
        exons=exons,
        protein_sequence=protein,
        cds_sequence=cds,
        protein_length=len(protein),
        exon_count=len(exons),
        completeness=completeness,
        hsps=hsps,
        reference_hit_count=len({hsp.query_id for hsp in hsps}),
    )
    candidate.tm_helix_count = predict_tm_helices(protein)

    if gff_path and gff_path.exists():
        overlaps = gff_features_overlapping(gff_path, seqid, locus_start, locus_end)
        if overlaps:
            candidate.validation_flags.append("existing_gff_overlap")

    if work_dir is not None:
        query_path = work_dir / f"{candidate.candidate_id}.fasta"
        write_fasta_records({candidate.candidate_id: protein}, query_path)
        ref_path = work_dir / TRUSTED_FASTA_NAME
        if ref_path.exists():
            identity, coverage, subject = mmseqs_best_hit(query_path, ref_path, work_dir)
            candidate.reference_identity = identity
            candidate.reference_coverage = coverage
            candidate.best_reference_hit = subject

    score_candidate(candidate)
    return candidate


def write_candidate_gff3(candidates: list[CandidateLocus], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    lines = ["##gff-version 3"]
    for cand in candidates:
        if cand.rank_category == RANK_REJECTED:
            continue
        gene_id = cand.candidate_id
        mrna_id = f"{gene_id}_mRNA"
        lines.append(
            "\t".join(
                [
                    cand.seqid,
                    f"{ACTIVE_FAMILY}_discovery",
                    "gene",
                    str(cand.locus_start),
                    str(cand.locus_end),
                    ".",
                    cand.strand,
                    ".",
                    f"ID={gene_id};Name={gene_id};rank={cand.rank_category};score={cand.rank_score}",
                ]
            )
        )
        lines.append(
            "\t".join(
                [
                    cand.seqid,
                    f"{ACTIVE_FAMILY}_discovery",
                    "mRNA",
                    str(cand.locus_start),
                    str(cand.locus_end),
                    ".",
                    cand.strand,
                    ".",
                    f"ID={mrna_id};Parent={gene_id}",
                ]
            )
        )
        for exon in cand.exons:
            exon_id = f"{gene_id}_exon{exon.exon_number}"
            lines.append(
                "\t".join(
                    [
                        cand.seqid,
                        f"{ACTIVE_FAMILY}_discovery",
                        "exon",
                        str(exon.start),
                        str(exon.end),
                        ".",
                        cand.strand,
                        ".",
                        f"ID={exon_id};Parent={mrna_id}",
                    ]
                )
            )
            lines.append(
                "\t".join(
                    [
                        cand.seqid,
                        f"{ACTIVE_FAMILY}_discovery",
                        "CDS",
                        str(exon.start),
                        str(exon.end),
                        ".",
                        cand.strand,
                        str(exon.phase),
                        f"ID={exon_id}_cds;Parent={mrna_id}",
                    ]
                )
            )
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_candidate_fasta(candidates: list[CandidateLocus], path: Path) -> None:
    records: dict[str, str] = {}
    for cand in candidates:
        if cand.rank_category == RANK_REJECTED:
            continue
        header = (
            f"{cand.candidate_id} "
            f"organism={cand.organism} "
            f"loc={cand.seqid}:{cand.locus_start}-{cand.locus_end} "
            f"rank={cand.rank_category} score={cand.rank_score} "
            f"tm={cand.tm_helix_count} exons={cand.exon_count}"
        )
        records[header] = cand.protein_sequence
    write_fasta_records(records, path)


def write_candidate_table(candidates: list[CandidateLocus], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "organism",
        "candidate_id",
        "seqid",
        "strand",
        "locus_start",
        "locus_end",
        "exon_count",
        "protein_length",
        "tm_helix_count",
        "reference_identity",
        "reference_coverage",
        "best_reference_hit",
        "reference_hit_count",
        "completeness",
        "rank_category",
        "rank_score",
        "rejection_reason",
        "validation_flags",
    ]
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for cand in candidates:
            writer.writerow(
                {
                    "organism": cand.organism,
                    "candidate_id": cand.candidate_id,
                    "seqid": cand.seqid,
                    "strand": cand.strand,
                    "locus_start": cand.locus_start,
                    "locus_end": cand.locus_end,
                    "exon_count": cand.exon_count,
                    "protein_length": cand.protein_length,
                    "tm_helix_count": cand.tm_helix_count,
                    "reference_identity": cand.reference_identity,
                    "reference_coverage": cand.reference_coverage,
                    "best_reference_hit": cand.best_reference_hit,
                    "reference_hit_count": cand.reference_hit_count,
                    "completeness": cand.completeness,
                    "rank_category": cand.rank_category,
                    "rank_score": cand.rank_score,
                    "rejection_reason": cand.rejection_reason,
                    "validation_flags": ";".join(cand.validation_flags),
                }
            )


def plot_candidate_exon_map(candidates: list[CandidateLocus], out_file: Path) -> None:
    import matplotlib.pyplot as plt
    from matplotlib.patches import Rectangle

    kept = [c for c in candidates if c.rank_category != RANK_REJECTED]
    if not kept:
        return

    fig_h = max(4, len(kept) * 0.55 + 1.5)
    fig, ax = plt.subplots(figsize=(14, fig_h))
    y = 0
    xmax = max(c.locus_end for c in kept)
    xmin = min(c.locus_start for c in kept)
    span = max(xmax - xmin, 1)

    for cand in kept:
        for exon in cand.exons:
            width = exon.end - exon.start + 1
            x = exon.start - xmin
            rect = Rectangle((x, y - 0.18), width, 0.36, facecolor="#9B59B6", edgecolor="black")
            ax.add_patch(rect)
            ax.text(
                x + width / 2,
                y,
                str(exon.exon_number),
                ha="center",
                va="center",
                fontsize=8,
                color="white",
                fontweight="bold",
            )
        label = f"{cand.candidate_id} ({cand.rank_category}, TM={cand.tm_helix_count})"
        ax.text(-span * 0.02, y, label, ha="right", va="center", fontsize=8)
        y += 1

    organism = kept[0].organism
    ax.set_title(f"Candidate {ACTIVE_FAMILY} exon map — {organism}")
    ax.set_xlabel("Genomic position")
    ax.set_yticks([])
    ax.set_xlim(-span * 0.35, span * 1.05)
    ax.set_ylim(-0.5, y)
    out_file.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(out_file, dpi=200)
    plt.close(fig)


def discover_for_species(
    organism: str,
    reference_fasta: Path,
    *,
    genome_dir: Path,
    results_dir: Path,
    skip_download: bool = False,
    force_redownload: bool = False,
) -> list[CandidateLocus]:
    species_results = results_dir / normalize_species_name(organism)
    species_results.mkdir(parents=True, exist_ok=True)
    work_dir = species_results / "_work"
    work_dir.mkdir(parents=True, exist_ok=True)

    genome_path = None
    if skip_download:
        local = sorted((genome_dir / normalize_species_name(organism)).glob("*genomic.fna*"))
        if not local:
            raise FileNotFoundError(f"No local genome for {organism} in {genome_dir}")
        genome_path = decompress_if_gzip(local[0])
    else:
        genome_path = ensure_genome_fasta(
            organism, genome_dir, force_redownload=force_redownload
        )

    miniprot_gff = work_dir / "miniprot.gff"
    run_miniprot_genome_search(genome_path, reference_fasta, miniprot_gff)
    raw_hsps = parse_miniprot_gff(miniprot_gff)
    hsps = filter_hsps(raw_hsps)
    loci = cluster_hsps_into_loci(hsps)

    gff_path = find_annotation_gff(organism)
    shutil.copy(reference_fasta, work_dir / TRUSTED_FASTA_NAME)

    candidates: list[CandidateLocus] = []
    for index, locus_hsps in enumerate(loci, start=1):
        candidate = build_candidate_locus(
            organism,
            index,
            locus_hsps,
            genome_path,
            gff_path=gff_path,
            work_dir=work_dir,
        )
        if candidate:
            candidates.append(candidate)

    candidates.sort(key=lambda item: item.rank_score, reverse=True)

    write_candidate_table(candidates, species_results / "candidate_loci.csv")
    write_candidate_gff3(candidates, species_results / CANDIDATE_GFF3_NAME)
    write_candidate_fasta(candidates, species_results / CANDIDATE_FASTA_NAME)
    try:
        plot_candidate_exon_map(
            candidates, species_results / "plots" / PLOT_NAME
        )
    except ImportError:
        pass

    summary = {
        "organism": organism,
        "genome_path": str(genome_path),
        "hsp_count_raw": len(raw_hsps),
        "hsp_count": len(hsps),
        "locus_count": len(loci),
        "candidate_count": len(candidates),
        "accepted_count": sum(1 for c in candidates if c.rank_category != RANK_REJECTED),
    }
    (species_results / "discovery_summary.json").write_text(
        json.dumps(summary, indent=2), encoding="utf-8"
    )
    return candidates


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Discover candidate gap-junction genes (innexin or connexin) in poorly "
            "annotated species using miniprot genome search, locus clustering, gene "
            "prediction, validation, and ranking."
        )
    )
    parser.add_argument(
        "--family",
        choices=sorted(FAMILY_DEFAULTS),
        default="innexin",
        help="Gap-junction gene family to discover (default: innexin).",
    )
    parser.add_argument(
        "--species-txt",
        type=Path,
        default=None,
        help="Species list text file (default: family-specific not_annotated_*.txt).",
    )
    parser.add_argument(
        "species",
        nargs="*",
        help="Optional species names (overrides / adds to species txt).",
    )
    parser.add_argument(
        "--reference-dir",
        type=Path,
        default=None,
        help="Directory with trusted family reference FASTAs.",
    )
    parser.add_argument(
        "--genome-dir",
        type=Path,
        default=None,
        help="Directory for downloaded target genomes.",
    )
    parser.add_argument(
        "--results-dir",
        type=Path,
        default=None,
        help="Output directory for discovery results.",
    )
    parser.add_argument(
        "--skip-download",
        action="store_true",
        help="Use existing genome FASTA only; do not download from NCBI.",
    )
    parser.add_argument(
        "--force-redownload",
        action="store_true",
        help="Delete existing genome files and download fresh copies from NCBI.",
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Process only the first N species.",
    )
    parser.add_argument(
        "--num-threads",
        type=int,
        default=max(1, os.cpu_count() or 1),
        help="CPU threads for miniprot/MMseqs2 (default: all cores).",
    )
    parser.add_argument(
        "--chunk-workers",
        type=int,
        default=CHUNKED_MAX_WORKERS,
        help="Parallel miniprot jobs for large genomes (default: 3).",
    )
    parser.add_argument(
        "--chunk-min-bp",
        type=int,
        default=CHUNKED_MIN_CONTIG_BP,
        help="Min contig size for chunked miniprot on large genomes.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    global SEARCH_THREADS, CHUNKED_MIN_CONTIG_BP, CHUNKED_MAX_WORKERS
    args = parse_args(argv)
    configure_family(args.family)
    if args.species_txt is None:
        args.species_txt = DEFAULT_SPECIES_TXT
    if args.reference_dir is None:
        args.reference_dir = REFERENCES_DIR
    if args.genome_dir is None:
        args.genome_dir = GENOME_DIR
    if args.results_dir is None:
        args.results_dir = RESULTS_DIR
    SEARCH_THREADS = max(1, args.num_threads)
    CHUNKED_MIN_CONTIG_BP = max(1, args.chunk_min_bp)
    CHUNKED_MAX_WORKERS = max(1, args.chunk_workers)
    require_discovery_tools()

    species = read_species_txt(args.species_txt) if args.species_txt.exists() else []
    species.extend(args.species)
    species = list(dict.fromkeys(species))
    if args.limit is not None:
        species = species[: args.limit]
    if not species:
        print("No species provided.")
        return 1

    reference_fasta, ref_meta = collect_reference_sequences(args.reference_dir)
    print(
        f"Loaded {len(ref_meta)} trusted {ACTIVE_FAMILY} references from {args.reference_dir}"
    )
    print(
        f"Using {SEARCH_THREADS} thread(s) for miniprot; "
        f"chunked mode: >= {CHUNKED_MIN_CONTIG_BP} bp, {CHUNKED_MAX_WORKERS} workers"
    )

    for index, organism in enumerate(species, start=1):
        species_dir = args.results_dir / normalize_species_name(organism)
        if (species_dir / "discovery_summary.json").exists():
            print(f"[{index}/{len(species)}] skip {organism} (already complete)")
            continue
        print(f"[{index}/{len(species)}] Discovering {ACTIVE_FAMILY}s in {organism}", flush=True)
        try:
            candidates = discover_for_species(
                organism,
                reference_fasta,
                genome_dir=args.genome_dir,
                results_dir=args.results_dir,
                skip_download=args.skip_download,
                force_redownload=args.force_redownload,
            )
            accepted = [c for c in candidates if c.rank_category != RANK_REJECTED]
            print(
                f"  -> {len(candidates)} candidates "
                f"({len(accepted)} accepted after validation/ranking)"
            )
        except Exception as exc:  # noqa: BLE001
            print(f"  warning: failed for {organism}: {exc}")

    print(f"Results saved under {args.results_dir}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
