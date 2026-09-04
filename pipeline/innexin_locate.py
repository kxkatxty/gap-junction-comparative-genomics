#!/usr/bin/env python3
"""Improved innexin locate: MMseqs2 region list → per-region miniprot.

miniprot alone returns roughly one best model per query protein, so it is a
floor on copy number. This module first finds all homologous genomic windows
with MMseqs2 (translated nucleotide search), then runs miniprot inside each
window and applies exon-count / length / chimera guards.
"""

from __future__ import annotations

import csv
import shutil
import subprocess
from dataclasses import asdict, dataclass
from pathlib import Path

from pipeline.batch_innexin_curator_probe import (
    _tool,
    classify_product,
    parse_miniprot_models,
    run_miniprot,
    translate_model,
)
from pipeline.common import PROJECT_ROOT

MMSEQS = "mmseqs"
TBLASTN = "tblastn"
MAKEBLASTDB = "makeblastdb"
DEFAULT_THREADS = 2
DEFAULT_EVALUE = 1e-5
DEFAULT_FLANK_BP = 5_000
DEFAULT_MERGE_GAP_BP = 2_000
# Tight NMS so neighbouring paralogs (ogre / Inx7, ~0.3 kb apart) stay separate;
# final 5 kb locus dedupe collapses redundant models of the same gene.
DEFAULT_NMS_SUPPRESS_BP = 250
DEFAULT_SEED_PAD_BP = 500
DEFAULT_MAX_EXONS = 20
DEFAULT_MAX_GENOMIC_SPAN_BP = 200_000
DEFAULT_MAX_SPAN_PER_AA = 400.0  # chimera guard: genomic span / protein aa


@dataclass
class Region:
    seqid: str
    start: int
    end: int
    best_evalue: float
    n_hits: int


@dataclass
class LocatedLocus:
    seqid: str
    start: int
    end: int
    strand: str
    identity: float
    target: str
    aa_length: int
    cys: int
    stops: int
    verdict: str
    exon_count: int
    genomic_span: int
    region_id: str
    method: str
    rejected_reason: str = ""


def ensure_fai(genome: Path) -> Path:
    fai = Path(str(genome) + ".fai")
    if not fai.exists():
        subprocess.run([_tool("samtools"), "faidx", str(genome)], check=True)
    return fai


def contig_lengths(genome: Path) -> dict[str, int]:
    ensure_fai(genome)
    lengths: dict[str, int] = {}
    with Path(str(genome) + ".fai").open() as fh:
        for line in fh:
            parts = line.split("\t")
            if len(parts) >= 2:
                lengths[parts[0]] = int(parts[1])
    return lengths


def _hits_from_blast_or_mmseqs_m8(m8: Path) -> list[tuple[str, int, int, float]]:
    """Parse BLAST/MMseqs tabular hits → (seqid, lo, hi, evalue)."""
    raw: list[tuple[str, int, int, float]] = []
    if not m8.exists():
        return raw
    for line in m8.read_text(encoding="utf-8", errors="replace").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 12:
            continue
        seqid = parts[1]
        sstart, send = int(parts[8]), int(parts[9])
        lo, hi = (sstart, send) if sstart <= send else (send, sstart)
        evalue_hit = float(parts[10])
        raw.append((seqid, lo, hi, evalue_hit))
    return raw


def run_tblastn_regions(
    genome: Path,
    queries: Path,
    work_dir: Path,
    *,
    threads: int = DEFAULT_THREADS,
    evalue: float = DEFAULT_EVALUE,
    merge_gap_bp: int = DEFAULT_MERGE_GAP_BP,
    nms_suppress_bp: int = DEFAULT_NMS_SUPPRESS_BP,
    seed_pad_bp: int = DEFAULT_SEED_PAD_BP,
) -> list[Region]:
    """Protein vs genome tblastn → NMS-seeded genomic regions (preferred on large genomes)."""
    work_dir.mkdir(parents=True, exist_ok=True)
    db = work_dir / "blast_db"
    m8 = work_dir / "tblastn_genome.m8"
    log = work_dir / "tblastn_genome.log"

    # Reuse BLAST DB if already built for this genome path.
    db_ok = (Path(str(db) + ".nsq").exists() or Path(str(db) + ".nal").exists()
             or Path(str(db) + ".nin").exists())
    if not db_ok:
        with log.open("w", encoding="utf-8") as log_handle:
            proc = subprocess.run(
                [
                    _tool(MAKEBLASTDB),
                    "-in",
                    str(genome),
                    "-dbtype",
                    "nucl",
                    "-out",
                    str(db),
                ],
                stdout=log_handle,
                stderr=subprocess.STDOUT,
            )
        if proc.returncode != 0:
            raise RuntimeError(f"makeblastdb failed (rc={proc.returncode}); see {log}")

    with log.open("a", encoding="utf-8") as log_handle:
        proc = subprocess.run(
            [
                _tool(TBLASTN),
                "-query",
                str(queries),
                "-db",
                str(db),
                "-evalue",
                str(evalue),
                "-outfmt",
                "6",
                "-num_threads",
                str(threads),
                "-max_target_seqs",
                "100",
                "-out",
                str(m8),
            ],
            stdout=log_handle,
            stderr=subprocess.STDOUT,
        )
    if proc.returncode != 0:
        raise RuntimeError(f"tblastn failed (rc={proc.returncode}); see {log}")

    raw = _hits_from_blast_or_mmseqs_m8(m8)
    seeds = seed_regions_nms(
        raw, suppress_bp=nms_suppress_bp, seed_pad_bp=seed_pad_bp
    )
    if seeds:
        return seeds
    return merge_hit_intervals(raw, merge_gap_bp=merge_gap_bp)


def run_mmseqs_regions(
    genome: Path,
    queries: Path,
    work_dir: Path,
    *,
    threads: int = DEFAULT_THREADS,
    evalue: float = DEFAULT_EVALUE,
    merge_gap_bp: int = DEFAULT_MERGE_GAP_BP,
    nms_suppress_bp: int = DEFAULT_NMS_SUPPRESS_BP,
    seed_pad_bp: int = DEFAULT_SEED_PAD_BP,
) -> list[Region]:
    """Protein vs genome translated MMseqs2 search → NMS-seeded regions.

    Note: on multi-GB genomes ORF extraction is slow; prefer run_tblastn_regions.
    """
    work_dir.mkdir(parents=True, exist_ok=True)
    m8 = work_dir / "mmseqs_genome.m8"
    tmp = work_dir / "mmseqs_tmp"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True, exist_ok=True)

    cmd = [
        _tool(MMSEQS),
        "easy-search",
        str(queries),
        str(genome),
        str(m8),
        str(tmp),
        "--search-type",
        "2",  # translated nucleotide target
        "-e",
        str(evalue),
        "--threads",
        str(threads),
        "--max-seqs",
        "100",
        "-s",
        "5.0",
    ]
    log = work_dir / "mmseqs_genome.log"
    with log.open("w", encoding="utf-8") as log_handle:
        proc = subprocess.run(cmd, stdout=log_handle, stderr=subprocess.STDOUT)
    if proc.returncode != 0:
        raise RuntimeError(f"mmseqs easy-search failed (rc={proc.returncode}); see {log}")

    raw = _hits_from_blast_or_mmseqs_m8(m8)
    seeds = seed_regions_nms(
        raw, suppress_bp=nms_suppress_bp, seed_pad_bp=seed_pad_bp
    )
    if seeds:
        return seeds
    return merge_hit_intervals(raw, merge_gap_bp=merge_gap_bp)


def seed_regions_nms(
    hits: list[tuple[str, int, int, float]],
    *,
    suppress_bp: int,
    seed_pad_bp: int,
) -> list[Region]:
    """Greedy non-maximum suppression on MMseqs HSPs (best e-value first)."""
    ordered = sorted(hits, key=lambda h: (h[3], -(h[2] - h[1])))
    kept: list[Region] = []
    for seqid, lo, hi, evalue in ordered:
        conflict = False
        for prev in kept:
            if prev.seqid != seqid:
                continue
            if not (hi + suppress_bp < prev.start or prev.end + suppress_bp < lo):
                conflict = True
                break
        if conflict:
            continue
        kept.append(
            Region(
                seqid=seqid,
                start=max(1, lo - seed_pad_bp),
                end=hi + seed_pad_bp,
                best_evalue=evalue,
                n_hits=1,
            )
        )
    kept.sort(key=lambda r: (r.seqid, r.start))
    return kept


def merge_hit_intervals(
    hits: list[tuple[str, int, int, float]],
    *,
    merge_gap_bp: int,
) -> list[Region]:
    by_contig: dict[str, list[tuple[int, int, float]]] = {}
    for seqid, lo, hi, evalue in hits:
        by_contig.setdefault(seqid, []).append((lo, hi, evalue))

    regions: list[Region] = []
    for seqid, intervals in by_contig.items():
        intervals.sort(key=lambda x: x[0])
        cur_lo, cur_hi, cur_e, n = intervals[0][0], intervals[0][1], intervals[0][2], 1
        for lo, hi, evalue in intervals[1:]:
            if lo <= cur_hi + merge_gap_bp:
                cur_hi = max(cur_hi, hi)
                cur_e = min(cur_e, evalue)
                n += 1
            else:
                regions.append(Region(seqid, cur_lo, cur_hi, cur_e, n))
                cur_lo, cur_hi, cur_e, n = lo, hi, evalue, 1
        regions.append(Region(seqid, cur_lo, cur_hi, cur_e, n))
    regions.sort(key=lambda r: (r.seqid, r.start))
    return regions


def extract_region_fasta(
    genome: Path,
    region: Region,
    out_fa: Path,
    *,
    flank_bp: int,
    lengths: dict[str, int],
) -> tuple[int, int]:
    """Write a flanked region FASTA. Returns (extract_start, extract_end) 1-based."""
    contig_len = lengths.get(region.seqid)
    if contig_len is None:
        raise KeyError(f"contig {region.seqid} missing from {genome}.fai")
    lo = max(1, region.start - flank_bp)
    hi = min(contig_len, region.end + flank_bp)
    ensure_fai(genome)
    proc = subprocess.run(
        [_tool("samtools"), "faidx", str(genome), f"{region.seqid}:{lo}-{hi}"],
        check=True,
        capture_output=True,
        text=True,
    )
    out_fa.parent.mkdir(parents=True, exist_ok=True)
    # Rewrite header to a simple id so miniprot coords stay local
    seq = "".join(
        line.strip() for line in proc.stdout.splitlines() if line and not line.startswith(">")
    )
    out_fa.write_text(f">{region.seqid}_{lo}_{hi}\n{seq}\n", encoding="utf-8")
    return lo, hi


def chimera_or_exon_reject(
    model: dict,
    prot: str,
    *,
    max_exons: int = DEFAULT_MAX_EXONS,
    max_span_bp: int = DEFAULT_MAX_GENOMIC_SPAN_BP,
    max_span_per_aa: float = DEFAULT_MAX_SPAN_PER_AA,
) -> str:
    exons = model.get("cds") or []
    n_exons = len(exons) if exons else 1
    span = int(model["end"]) - int(model["start"]) + 1
    clean = prot.split("*")[0]
    aa = max(len(clean), 1)
    if n_exons > max_exons:
        return f"too_many_exons:{n_exons}"
    if span > max_span_bp:
        return f"genomic_span_too_large:{span}"
    if span / aa > max_span_per_aa:
        return f"chimera_span_per_aa:{span / aa:.1f}"
    return ""


def locate_with_mmseqs_then_miniprot(
    genome: Path,
    queries: Path,
    out_dir: Path,
    *,
    threads: int = DEFAULT_THREADS,
    evalue: float = DEFAULT_EVALUE,
    flank_bp: int = DEFAULT_FLANK_BP,
    extra_flanks: tuple[int, ...] = (30_000,),
    merge_gap_bp: int = DEFAULT_MERGE_GAP_BP,
    outs: float = 0.5,
    max_exons: int = DEFAULT_MAX_EXONS,
    region_finder: str = "tblastn",
) -> list[LocatedLocus]:
    """Region list → per-region miniprot locate pipeline.

    region_finder:
      - "tblastn" (default): sturdy + fast on large genomes
      - "mmseqs": translated MMseqs2 (fine for small genomes; slow ORF extract on multi-GB)
    """
    out_dir.mkdir(parents=True, exist_ok=True)
    lengths = contig_lengths(genome)
    finder = region_finder.lower().strip()
    if finder == "tblastn":
        regions = run_tblastn_regions(
            genome,
            queries,
            out_dir / "tblastn_work",
            threads=threads,
            evalue=evalue,
            merge_gap_bp=merge_gap_bp,
        )
        method_tag = "tblastn-regions+miniprot"
    elif finder in {"mmseqs", "mmseqs2"}:
        regions = run_mmseqs_regions(
            genome,
            queries,
            out_dir / "mmseqs_work",
            threads=threads,
            evalue=evalue,
            merge_gap_bp=merge_gap_bp,
        )
        method_tag = "mmseqs-regions+miniprot"
    else:
        raise ValueError(f"Unknown region_finder={region_finder!r}")
    _write_regions_csv(out_dir / "regions.csv", regions)

    loci: list[LocatedLocus] = []
    regions_dir = out_dir / "regions"
    if regions_dir.exists():
        shutil.rmtree(regions_dir)
    regions_dir.mkdir(parents=True, exist_ok=True)

    flanks = []
    for f in (flank_bp, *extra_flanks):
        if f not in flanks:
            flanks.append(f)

    for i, region in enumerate(regions):
        for flank in flanks:
            rid = f"R{i:04d}_f{flank}"
            region_fa = regions_dir / f"{rid}.fna"
            gff = regions_dir / f"{rid}.gff"
            log = regions_dir / f"{rid}.miniprot.log"
            extract_start, _extract_end = extract_region_fasta(
                genome, region, region_fa, flank_bp=flank, lengths=lengths
            )
            rc = run_miniprot(region_fa, queries, gff, log, threads, outs=outs)
            if rc != 0:
                continue
            models = parse_miniprot_models(gff)
            for model in models:
                local_start = int(model["start"])
                local_end = int(model["end"])
                genome_start = extract_start + local_start - 1
                genome_end = extract_start + local_end - 1
                genome_model = {
                    **model,
                    "seqid": region.seqid,
                    "start": genome_start,
                    "end": genome_end,
                    "cds": [
                        (extract_start + a - 1, extract_start + b - 1)
                        for a, b in (model.get("cds") or [])
                    ],
                }
                try:
                    prot = translate_model(genome, genome_model)
                except (subprocess.CalledProcessError, ValueError, KeyError):
                    continue
                reason = chimera_or_exon_reject(
                    genome_model, prot, max_exons=max_exons
                )
                verdict, motif = classify_product(float(genome_model["identity"]), prot)
                if reason:
                    verdict = "rejected"
                elif verdict == "noise":
                    continue
                trimmed = prot[:-1] if prot.endswith("*") else prot
                clean = trimmed.split("*")[0]
                loci.append(
                    LocatedLocus(
                        seqid=region.seqid,
                        start=genome_start,
                        end=genome_end,
                        strand=str(genome_model["strand"]),
                        identity=float(genome_model["identity"]),
                        target=str(genome_model.get("target") or ""),
                        aa_length=len(clean),
                        cys=clean.count("C"),
                        stops=trimmed.count("*"),
                        verdict=verdict,
                        exon_count=len(genome_model.get("cds") or []) or 1,
                        genomic_span=genome_end - genome_start + 1,
                        region_id=rid,
                        method=method_tag,
                        rejected_reason=reason or motif.detail,
                    )
                )

    # Prefer present > fragmentary > weak, then longer aa, then higher identity
    kept = dedupe_loci(loci, window_bp=5000)
    _write_loci_csv(out_dir / "loci.csv", kept)
    return kept


def locate_miniprot_only(
    genome: Path,
    queries: Path,
    out_dir: Path,
    *,
    threads: int = DEFAULT_THREADS,
    outs: float = 0.5,
) -> list[LocatedLocus]:
    """Current curator-style locate (floor): whole-genome miniprot + 5 kb dedupe."""
    out_dir.mkdir(parents=True, exist_ok=True)
    gff = out_dir / "miniprot.gff"
    log = out_dir / "miniprot.log"
    rc = run_miniprot(genome, queries, gff, log, threads, outs=outs)
    if rc != 0:
        raise RuntimeError(f"miniprot failed rc={rc}; see {log}")
    models = parse_miniprot_models(gff)
    best: dict[tuple[str, int], dict] = {}
    for model in models:
        key = (model["seqid"], model["start"] // 5000)
        prev = best.get(key)
        if prev is None or model["identity"] > prev["identity"]:
            best[key] = model

    loci: list[LocatedLocus] = []
    for model in best.values():
        try:
            prot = translate_model(genome, model)
        except (subprocess.CalledProcessError, ValueError, KeyError):
            continue
        verdict, motif = classify_product(model["identity"], prot)
        if verdict == "noise":
            continue
        clean = prot.split("*")[0]
        loci.append(
            LocatedLocus(
                seqid=model["seqid"],
                start=int(model["start"]),
                end=int(model["end"]),
                strand=str(model["strand"]),
                identity=float(model["identity"]),
                target=str(model.get("target") or ""),
                aa_length=len(clean),
                cys=clean.count("C"),
                stops=prot.count("*"),
                verdict=verdict,
                exon_count=len(model.get("cds") or []) or 1,
                genomic_span=int(model["end"]) - int(model["start"]) + 1,
                region_id="genome",
                method="miniprot_only",
                rejected_reason=motif.detail,
            )
        )
    _write_loci_csv(out_dir / "loci.csv", loci)
    return loci


def dedupe_loci(loci: list[LocatedLocus], *, window_bp: int = 5000) -> list[LocatedLocus]:
    rank = {"present": 3, "fragmentary": 2, "weak": 1, "rejected": 0}

    def better(a: LocatedLocus, b: LocatedLocus) -> bool:
        """Return True if a should replace b."""
        ra, rb = rank.get(a.verdict, 0), rank.get(b.verdict, 0)
        if ra != rb:
            return ra > rb
        if a.aa_length != b.aa_length:
            return a.aa_length > b.aa_length
        return a.identity > b.identity

    best: dict[tuple, LocatedLocus] = {}
    for loc in loci:
        if loc.verdict == "rejected":
            key = (loc.seqid, loc.start // window_bp, "rejected", loc.start)
            prev = best.get(key)
            if prev is None or better(loc, prev):
                best[key] = loc
            continue
        key = (loc.seqid, loc.start // window_bp)
        prev = best.get(key)
        if prev is None or better(loc, prev):
            best[key] = loc
    return sorted(best.values(), key=lambda x: (x.seqid, x.start))


def present_loci(loci: list[LocatedLocus]) -> list[LocatedLocus]:
    return [x for x in loci if x.verdict == "present"]


def _write_regions_csv(path: Path, regions: list[Region]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(
            fh, fieldnames=["seqid", "start", "end", "best_evalue", "n_hits"]
        )
        w.writeheader()
        for r in regions:
            w.writerow(asdict(r))


def _write_loci_csv(path: Path, loci: list[LocatedLocus]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not loci:
        path.write_text(
            "seqid,start,end,strand,identity,target,aa_length,cys,stops,verdict,"
            "exon_count,genomic_span,region_id,method,rejected_reason\n",
            encoding="utf-8",
        )
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(asdict(loci[0]).keys()))
        w.writeheader()
        for loc in loci:
            row = asdict(loc)
            row["identity"] = f"{loc.identity:.4f}"
            w.writerow(row)


def default_curator_queries() -> Path:
    from pipeline.batch_innexin_curator_probe import QUERY_FASTA, build_query_fasta

    if QUERY_FASTA.exists() and QUERY_FASTA.stat().st_size > 0:
        return QUERY_FASTA
    return build_query_fasta(QUERY_FASTA)


def dmel_reference_queries() -> Path:
    """All 8 Dmel innexin proteins (for optional self-query tests)."""
    ref_dir = (
        PROJECT_ROOT
        / "project/data/references/innexins/Drosophila_melanogaster"
    )
    out = (
        PROJECT_ROOT
        / "project/results/dmel_innexin_benchmark/queries_dmel8.fasta"
    )
    out.parent.mkdir(parents=True, exist_ok=True)
    chunks: list[str] = []
    for path in sorted(ref_dir.glob("*.fasta")):
        chunks.append(path.read_text(encoding="utf-8"))
    out.write_text("".join(chunks), encoding="utf-8")
    return out
