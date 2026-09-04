#!/usr/bin/env python3
"""Batch genome-direct innexin probe (gene-curator style) across discovery genomes.

For each species genome:
  1. miniprot trusted innexin queries (coarse locate)
  2. translate CDS products
  3. classify: present / fragmentary / weak / none

Writes a searched/found summary table at the end.
"""

from __future__ import annotations

import argparse
import csv
import json
import subprocess
import sys
from dataclasses import dataclass
from pathlib import Path

from Bio.Seq import Seq

from pipeline.common import PROJECT_ROOT, RESULTS_DIR, write_csv

GENOME_ROOT = PROJECT_ROOT / "project" / "data" / "genomes" / "innexin_discovery"
DISCOVERY_ROOT = RESULTS_DIR / "innexin_discovery"
OUT_ROOT = RESULTS_DIR / "gene_curator_probe"
QUERY_FASTA = OUT_ROOT / "batch_queries.fasta"
QUERY_FASTA_EXPANDED = OUT_ROOT / "batch_queries_expanded.fasta"
SUMMARY_CSV = OUT_ROOT / "innexin_search_summary.csv"
SUMMARY_MD = OUT_ROOT / "innexin_search_summary.md"

MINIPROT = "miniprot"
SAMTOOLS = "samtools"
SYNVOY_BIN = Path.home() / "miniconda3" / "envs" / "synvoy_env" / "bin"


@dataclass
class Hit:
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
    el1_cys: str = ""
    el2_cys: str = ""
    el_motif: str = ""  # pass|fail|unresolved
    tm_pred: str = ""


def _tool(name: str) -> str:
    cand = SYNVOY_BIN / name
    return str(cand) if cand.exists() else name


def build_query_fasta(path: Path) -> Path:
    path.parent.mkdir(parents=True, exist_ok=True)
    chunks: list[str] = []
    # Rotifer positive-control products if present
    adineta = DISCOVERY_ROOT / "Adineta_vaga" / "candidate_innexins.fasta"
    if adineta.exists() and adineta.stat().st_size > 0:
        chunks.append(adineta.read_text(encoding="utf-8"))
    # Full Dmel set (production): covers paralogs that a 3-query pack misses on the
    # Drosophila ground-truth benchmark (Inx3/Inx5/Inx6/Inx7/zpg).
    refs = [
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/ogre__P27716.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/Inx2__Q9V427.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/Inx3__Q9VAS7.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/zpg__Q9VRX6.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/Inx5__Q9VWL5.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/Inx6__Q9VR82.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/Inx7__Q9V3W6.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/shakB__P33085.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Caenorhabditis_elegans/inx-2__Q9U3K5.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Caenorhabditis_elegans/unc-9__O01393.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Caenorhabditis_elegans/unc-7__Q03412.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Aedes_aegypti/shakB__Q1DH70.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Schistocerca_americana/inx2__Q9XYN1.fasta",
    ]
    for ref in refs:
        if ref.exists():
            chunks.append(ref.read_text(encoding="utf-8"))
    if not chunks:
        raise FileNotFoundError("No innexin query FASTA sources found")
    path.write_text("".join(chunks), encoding="utf-8")
    return path


def find_genome_fasta(species_dir: Path) -> Path | None:
    fnas = list(species_dir.glob("*.fna")) + list(species_dir.glob("*.fa")) + list(species_dir.glob("*.fasta"))
    if not fnas:
        return None
    decompressed = [p for p in fnas if "decompressed" in p.name]
    pool = decompressed or fnas
    return max(pool, key=lambda p: p.stat().st_size)


def prior_discovery(species: str) -> tuple[str, str]:
    summary = DISCOVERY_ROOT / species / "discovery_summary.json"
    if not summary.exists():
        return "not_run", ""
    data = json.loads(summary.read_text(encoding="utf-8"))
    return (
        f"accepted={data.get('accepted_count', '?')};candidates={data.get('candidate_count', '?')}",
        str(data.get("accepted_count", "")),
    )


def ensure_fai(genome: Path) -> None:
    fai = Path(str(genome) + ".fai")
    if fai.exists():
        return
    subprocess.run([_tool(SAMTOOLS), "faidx", str(genome)], check=True, capture_output=True)


def run_miniprot(
    genome: Path,
    queries: Path,
    out_gff: Path,
    log_path: Path,
    threads: int,
    *,
    outs: float = 0.5,
) -> int:
    out_gff.parent.mkdir(parents=True, exist_ok=True)
    cmd = [
        _tool(MINIPROT),
        "-t",
        str(threads),
        f"--outs={outs}",
        "--gff",
        str(genome),
        str(queries),
    ]
    with out_gff.open("w", encoding="utf-8") as gff_handle, log_path.open(
        "w", encoding="utf-8"
    ) as log_handle:
        proc = subprocess.run(cmd, stdout=gff_handle, stderr=log_handle)
    return proc.returncode


def parse_miniprot_models(gff: Path) -> list[dict]:
    models: dict[str, dict] = {}
    current = None
    if not gff.exists():
        return []
    for line in gff.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith("#") or "\t" not in line:
            continue
        parts = line.split("\t")
        if len(parts) < 9:
            continue
        attrs = dict(x.split("=", 1) for x in parts[8].split(";") if "=" in x)
        if parts[2] == "mRNA":
            mid = attrs.get("ID", f"MP{len(models)}")
            current = mid
            models[mid] = {
                "id": mid,
                "seqid": parts[0],
                "start": int(parts[3]),
                "end": int(parts[4]),
                "strand": parts[6],
                "identity": float(attrs.get("Identity") or 0),
                "target": attrs.get("Target", ""),
                "cds": [],
            }
        elif parts[2] == "CDS" and current and current in models:
            parent = attrs.get("Parent", current)
            if parent in models:
                models[parent]["cds"].append((int(parts[3]), int(parts[4])))
    return list(models.values())


def faidx_region(genome: Path, seqid: str, start: int, end: int) -> str:
    region = f"{seqid}:{start}-{end}"
    proc = subprocess.run(
        [_tool(SAMTOOLS), "faidx", str(genome), region],
        check=True,
        capture_output=True,
        text=True,
    )
    lines = proc.stdout.splitlines()
    return "".join(line.strip() for line in lines if line and not line.startswith(">"))


def translate_model(genome: Path, model: dict) -> str:
    cds = model["cds"]
    if not cds:
        return ""
    strand = model["strand"]
    ordered = sorted(cds, key=lambda x: x[0], reverse=(strand == "-"))
    pieces: list[str] = []
    for start, end in ordered:
        dna = faidx_region(genome, model["seqid"], start, end)
        if strand == "-":
            dna = str(Seq(dna).reverse_complement())
        pieces.append(dna)
    joined = "".join(pieces)
    if len(joined) % 3:
        joined = joined[: len(joined) - (len(joined) % 3)]
    if not joined:
        return ""
    return str(Seq(joined).translate(to_stop=False))


def classify_product(identity: float, prot: str) -> tuple[str, object]:
    """Return (verdict, ElCysMotif). Present prefers EL 2+2 over whole-protein Cys≥4."""
    from pipeline.innexin_topology import evaluate_el_cys_motif

    # Ignore a single terminal stop codon; count only internal stops.
    trimmed = prot[:-1] if prot.endswith("*") else prot
    clean = trimmed.split("*")[0]
    stops = trimmed.count("*")
    cys = clean.count("C")
    clean_aa = len(clean)
    motif = evaluate_el_cys_motif(clean)

    if identity >= 0.20 and clean_aa >= 280 and stops == 0:
        if motif.status == "pass":
            return "present", motif
        if motif.status == "fail":
            return "fragmentary", motif
        if cys >= 4:
            return "present", motif
        return "fragmentary", motif
    if identity >= 0.20 and clean_aa >= 180 and stops == 0:
        return "fragmentary", motif
    if identity >= 0.20 and clean_aa >= 100:
        return "weak", motif
    return "noise", motif


def assess_species(genome: Path, gff: Path) -> list[Hit]:
    models = parse_miniprot_models(gff)
    # Deduplicate overlapping models: keep highest identity per (seqid, rounded window)
    best: dict[tuple[str, int], dict] = {}
    for model in models:
        key = (model["seqid"], model["start"] // 5000)
        prev = best.get(key)
        if prev is None or model["identity"] > prev["identity"]:
            best[key] = model

    hits: list[Hit] = []
    for model in sorted(best.values(), key=lambda m: -m["identity"]):
        try:
            prot = translate_model(genome, model)
        except (subprocess.CalledProcessError, ValueError, KeyError):
            continue
        clean = prot.split("*")[0]
        verdict, motif = classify_product(model["identity"], prot)
        if verdict == "noise":
            continue
        hits.append(
            Hit(
                seqid=model["seqid"],
                start=model["start"],
                end=model["end"],
                strand=model["strand"],
                identity=model["identity"],
                target=model["target"],
                aa_length=len(clean),
                cys=clean.count("C"),
                stops=prot.count("*"),
                verdict=verdict,
                el1_cys="" if motif.el1_cys is None else str(motif.el1_cys),
                el2_cys="" if motif.el2_cys is None else str(motif.el2_cys),
                el_motif=motif.status,
                tm_pred=str(motif.tm_count),
            )
        )
    return hits


def iter_species(
    *,
    max_mb: float,
    include: set[str] | None,
    limit: int | None,
) -> list[tuple[str, Path]]:
    rows: list[tuple[str, Path, int]] = []
    for species_dir in sorted(GENOME_ROOT.iterdir()):
        if not species_dir.is_dir():
            continue
        if include and species_dir.name not in include:
            continue
        genome = find_genome_fasta(species_dir)
        if genome is None:
            continue
        size = genome.stat().st_size
        if size < 1_000_000:  # empty / stub
            continue
        if size > max_mb * 1_000_000:
            continue
        rows.append((species_dir.name, genome, size))
    rows.sort(key=lambda r: r[2])  # small genomes first
    if limit is not None:
        rows = rows[:limit]
    return [(name, path) for name, path, _ in rows]


def write_species_hits(out_dir: Path, hits: list[Hit]) -> None:
    rows = [
        {
            "seqid": h.seqid,
            "start": str(h.start),
            "end": str(h.end),
            "strand": h.strand,
            "identity": f"{h.identity:.4f}",
            "aa_length": str(h.aa_length),
            "cys": str(h.cys),
            "stops": str(h.stops),
            "verdict": h.verdict,
            "target": h.target,
            "el1_cys": h.el1_cys,
            "el2_cys": h.el2_cys,
            "el_motif": h.el_motif,
            "tm_pred": h.tm_pred,
        }
        for h in hits
    ]
    write_csv(out_dir / "hits.csv", rows)


def summarize(rows: list[dict[str, str]]) -> None:
    for row in rows:
        row.setdefault("locate_method", row.get("locate_method", ""))
        row.setdefault("notes", row.get("notes", ""))
    fieldnames: list[str] = []
    for row in rows:
        for key in row.keys():
            if key not in fieldnames:
                fieldnames.append(key)
    write_csv(SUMMARY_CSV, rows, fieldnames=fieldnames)
    present = [r for r in rows if int(r.get("present_loci") or 0) > 0]
    frag = [r for r in rows if int(r.get("fragmentary_loci") or 0) > 0 and int(r.get("present_loci") or 0) == 0]
    none = [r for r in rows if r.get("status") == "searched" and int(r.get("present_loci") or 0) == 0 and int(r.get("fragmentary_loci") or 0) == 0]
    failed = [r for r in rows if r.get("status") != "searched"]
    total_present = sum(int(r.get("present_loci") or 0) for r in present)
    lines = [
        "# Innexin curator batch search — summary",
        "",
        f"Species searched: **{sum(1 for r in rows if r.get('status')=='searched')}**",
        f"With ≥1 **present** locus: **{len(present)}**",
        f"Total present loci (**lower bound / floor**): **≥ {total_present}**",
        f"Only fragmentary/weak: **{len(frag)}**",
        f"No usable hit: **{len(none)}**",
        f"Failed / skipped: **{len(failed)}**",
        "",
        "Counts are **floors**, not full paralog censuses (Dmel calibration: "
        "whole-genome miniprot recovers 3/8 known innexins; region seeding "
        "+ per-region miniprot recovers 7–8/8 depending on query pack).",
        "",
        "## Present (likely innexin loci recovered; counts are ≥ floors)",
        "",
    ]
    for r in present:
        lines.append(
            f"- **{r['species']}**: present=≥{r['present_loci']}, fragmentary={r['fragmentary_loci']}, "
            f"weak={r['weak_loci']}, best_id={r['best_identity']}, best_aa={r['best_aa']} "
            f"(prior discovery: {r['prior_discovery']}; locate={r.get('locate_method','')})"
        )
    lines.extend(["", "## Fragmentary / weak only", ""])
    for r in frag:
        lines.append(
            f"- **{r['species']}**: fragmentary={r['fragmentary_loci']}, weak={r['weak_loci']}, "
            f"best_id={r['best_identity']}, best_aa={r['best_aa']}"
        )
    lines.extend(["", "## No usable hit", ""])
    for r in none:
        lines.append(f"- **{r['species']}** (miniprot mRNAs={r['miniprot_mrna_count']})")
    lines.extend(["", "## Failed / skipped", ""])
    for r in failed:
        lines.append(f"- **{r['species']}**: {r['status']} — {r.get('notes','')}")
    lines.extend(
        [
            "",
            f"Machine-readable table: `{SUMMARY_CSV.relative_to(PROJECT_ROOT)}`",
            "",
        ]
    )
    SUMMARY_MD.write_text("\n".join(lines) + "\n", encoding="utf-8")


def located_to_hits(loci) -> list[Hit]:
    """Convert innexin_locate.LocatedLocus rows into curator Hit objects."""
    hits: list[Hit] = []
    for loc in loci:
        if loc.verdict == "rejected":
            continue
        hits.append(
            Hit(
                seqid=loc.seqid,
                start=int(loc.start),
                end=int(loc.end),
                strand=str(loc.strand),
                identity=float(loc.identity),
                target=str(loc.target),
                aa_length=int(loc.aa_length),
                cys=int(loc.cys),
                stops=int(loc.stops),
                verdict=str(loc.verdict),
            )
        )
    return hits


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Batch innexin curator probe across discovery genomes")
    p.add_argument("--max-mb", type=float, default=2000.0, help="Skip genomes larger than this (MB)")
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--limit", type=int, default=None, help="Max species to process")
    p.add_argument("--species", action="append", default=[], help="Only these folder names")
    p.add_argument(
        "--species-file",
        type=Path,
        default=None,
        help="Text file with one species folder name per line",
    )
    p.add_argument("--resume", action="store_true", help="Skip species that already have hits.csv")
    p.add_argument("--force", action="store_true", help="Re-run even if hits.csv already exists")
    p.add_argument("--queries", type=Path, default=None, help="Query FASTA (default: build/use expanded if present)")
    p.add_argument("--outs", type=float, default=0.5, help="miniprot --outs threshold (lower = more sensitive)")
    p.add_argument("--use-expanded", action="store_true", help="Prefer batch_queries_expanded.fasta if present")
    p.add_argument(
        "--locate",
        choices=("tblastn-regions", "mmseqs-regions", "miniprot"),
        default="tblastn-regions",
        help="Locate method: tblastn region list → per-region miniprot (default), "
        "MMseqs2 regions → miniprot, or legacy whole-genome miniprot",
    )
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    OUT_ROOT.mkdir(parents=True, exist_ok=True)
    if args.queries is not None:
        queries = args.queries
        if not queries.exists():
            raise FileNotFoundError(queries)
    elif args.use_expanded and QUERY_FASTA_EXPANDED.exists():
        queries = QUERY_FASTA_EXPANDED
    else:
        # Always rebuild the production pack (full Dmel set + Cele/insect anchors).
        queries = build_query_fasta(QUERY_FASTA)
    include: set[str] | None = None
    names: list[str] = []
    if args.species:
        names.extend(args.species)
    if args.species_file is not None:
        names.extend(
            line.strip()
            for line in args.species_file.read_text(encoding="utf-8").splitlines()
            if line.strip() and not line.strip().startswith("#")
        )
    if names:
        include = {s.replace(" ", "_") for s in names}
    species_list = iter_species(max_mb=args.max_mb, include=include, limit=args.limit)
    if not species_list:
        print("No genomes matched filters.", file=sys.stderr)
        return 1

    print(f"Query pack: {queries} ({queries.stat().st_size} bytes)")
    print(
        f"Species to search: {len(species_list)} "
        f"(max {args.max_mb} MB, outs={args.outs}, locate={args.locate})"
    )

    summary_rows: list[dict[str, str]] = []
    # Keep previously completed rows when resume + rewriting summary
    existing: dict[str, dict[str, str]] = {}
    if SUMMARY_CSV.exists():
        with SUMMARY_CSV.open(encoding="utf-8", newline="") as handle:
            for row in csv.DictReader(handle):
                existing[row["species"]] = row

    for i, (species, genome) in enumerate(species_list, 1):
        out_dir = OUT_ROOT / species
        hits_path = out_dir / "hits.csv"
        prior, prior_accepted = prior_discovery(species)
        print(f"[{i}/{len(species_list)}] {species} ({genome.stat().st_size/1e6:.0f} MB) prior={prior}")

        if args.resume and not args.force and hits_path.exists():
            print("  SKIP resume")
            if species in existing:
                summary_rows.append(existing[species])
            continue

        out_dir.mkdir(parents=True, exist_ok=True)
        try:
            ensure_fai(genome)
            if args.locate in {"tblastn-regions", "mmseqs-regions"}:
                from pipeline.innexin_locate import locate_with_mmseqs_then_miniprot

                finder = "tblastn" if args.locate.startswith("tblastn") else "mmseqs"
                located = locate_with_mmseqs_then_miniprot(
                    genome,
                    queries,
                    out_dir / f"locate_{finder}_regions",
                    threads=args.threads,
                    outs=args.outs,
                    region_finder=finder,
                )
                hits = located_to_hits(located)
                n_models = len(located)
                locate_note = f"{finder}-regions+miniprot"
            else:
                gff = out_dir / "miniprot.gff"
                log = out_dir / "miniprot.log"
                rc = run_miniprot(genome, queries, gff, log, args.threads, outs=args.outs)
                if rc != 0:
                    summary_rows.append(
                        {
                            "species": species,
                            "status": "miniprot_failed",
                            "genome": str(genome.relative_to(PROJECT_ROOT)),
                            "prior_discovery": prior,
                            "prior_accepted": prior_accepted,
                            "miniprot_mrna_count": "0",
                            "present_loci": "0",
                            "fragmentary_loci": "0",
                            "weak_loci": "0",
                            "best_identity": "",
                            "best_aa": "",
                            "best_locus": "",
                            "locate_method": "miniprot",
                            "notes": f"exit={rc}",
                        }
                    )
                    continue
                n_models = len(parse_miniprot_models(gff))
                hits = assess_species(genome, gff)
                locate_note = "miniprot"

            write_species_hits(out_dir, hits)
            present = [h for h in hits if h.verdict == "present"]
            frag = [h for h in hits if h.verdict == "fragmentary"]
            weak = [h for h in hits if h.verdict == "weak"]
            best = max(hits, key=lambda h: (h.aa_length, h.identity), default=None)
            row = {
                "species": species,
                "status": "searched",
                "genome": str(genome.relative_to(PROJECT_ROOT)),
                "prior_discovery": prior,
                "prior_accepted": prior_accepted,
                "miniprot_mrna_count": str(n_models),
                "present_loci": str(len(present)),
                "fragmentary_loci": str(len(frag)),
                "weak_loci": str(len(weak)),
                "best_identity": f"{best.identity:.4f}" if best else "",
                "best_aa": str(best.aa_length) if best else "",
                "best_locus": f"{best.seqid}:{best.start}-{best.end}" if best else "",
                "locate_method": locate_note,
                "notes": "present_loci is a lower bound (floor)",
            }
            summary_rows.append(row)
            print(
                f"  models={n_models} present=≥{len(present)} frag={len(frag)} weak={len(weak)}"
                + (f" best={best.aa_length}aa id={best.identity:.2f}" if best else "")
                + f" [{locate_note}]"
            )
        except Exception as exc:  # noqa: BLE001 — batch must continue
            summary_rows.append(
                {
                    "species": species,
                    "status": "error",
                    "genome": str(genome.relative_to(PROJECT_ROOT)),
                    "prior_discovery": prior,
                    "prior_accepted": prior_accepted,
                    "miniprot_mrna_count": "0",
                    "present_loci": "0",
                    "fragmentary_loci": "0",
                    "weak_loci": "0",
                    "best_identity": "",
                    "best_aa": "",
                    "best_locus": "",
                    "locate_method": args.locate,
                    "notes": str(exc)[:200],
                }
            )
            print(f"  ERROR {exc}")

    # Merge any existing species not re-run
    seen = {r["species"] for r in summary_rows}
    for sp, row in existing.items():
        if sp not in seen:
            summary_rows.append(row)
    summary_rows.sort(key=lambda r: r["species"])
    summarize(summary_rows)
    print(f"\nSummary CSV: {SUMMARY_CSV}")
    print(f"Summary MD:  {SUMMARY_MD}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
