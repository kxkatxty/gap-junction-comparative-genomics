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
    refs = [
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/Inx2__Q9V427.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/shakB__P33085.fasta",
        PROJECT_ROOT / "project/data/references/innexins/Drosophila_melanogaster/ogre__P27716.fasta",
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


def classify_product(identity: float, prot: str) -> str:
    aa = len(prot.split("*")[0])
    stops = prot.count("*")
    cys = prot.split("*")[0].count("C")
    clean_aa = aa
    if identity >= 0.20 and clean_aa >= 280 and stops == 0 and cys >= 4:
        return "present"
    if identity >= 0.20 and clean_aa >= 180 and stops == 0:
        return "fragmentary"
    if identity >= 0.20 and clean_aa >= 100:
        return "weak"
    return "noise"


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
        verdict = classify_product(model["identity"], prot)
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
        }
        for h in hits
    ]
    write_csv(out_dir / "hits.csv", rows)
    fasta = out_dir / "present_or_fragmentary.fasta"
    # products are not stored here (would need re-translate); keep table only


def summarize(rows: list[dict[str, str]]) -> None:
    write_csv(SUMMARY_CSV, rows)
    present = [r for r in rows if int(r.get("present_loci") or 0) > 0]
    frag = [r for r in rows if int(r.get("fragmentary_loci") or 0) > 0 and int(r.get("present_loci") or 0) == 0]
    none = [r for r in rows if r.get("status") == "searched" and int(r.get("present_loci") or 0) == 0 and int(r.get("fragmentary_loci") or 0) == 0]
    failed = [r for r in rows if r.get("status") != "searched"]
    lines = [
        "# Innexin curator batch search — summary",
        "",
        f"Species searched: **{sum(1 for r in rows if r.get('status')=='searched')}**",
        f"With ≥1 **present** locus: **{len(present)}**",
        f"Only fragmentary/weak: **{len(frag)}**",
        f"No usable hit: **{len(none)}**",
        f"Failed / skipped: **{len(failed)}**",
        "",
        "## Present (likely innexin loci recovered)",
        "",
    ]
    for r in present:
        lines.append(
            f"- **{r['species']}**: present={r['present_loci']}, fragmentary={r['fragmentary_loci']}, "
            f"weak={r['weak_loci']}, best_id={r['best_identity']}, best_aa={r['best_aa']} "
            f"(prior discovery: {r['prior_discovery']})"
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
        queries = build_query_fasta(QUERY_FASTA)
        if QUERY_FASTA_EXPANDED.exists() and QUERY_FASTA_EXPANDED.stat().st_size > queries.stat().st_size:
            queries = QUERY_FASTA_EXPANDED
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
    print(f"Species to search: {len(species_list)} (max {args.max_mb} MB, outs={args.outs})")

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
        gff = out_dir / "miniprot.gff"
        log = out_dir / "miniprot.log"
        try:
            ensure_fai(genome)
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
                        "notes": f"exit={rc}",
                    }
                )
                continue
            models = parse_miniprot_models(gff)
            hits = assess_species(genome, gff)
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
                "miniprot_mrna_count": str(len(models)),
                "present_loci": str(len(present)),
                "fragmentary_loci": str(len(frag)),
                "weak_loci": str(len(weak)),
                "best_identity": f"{best.identity:.4f}" if best else "",
                "best_aa": str(best.aa_length) if best else "",
                "best_locus": f"{best.seqid}:{best.start}-{best.end}" if best else "",
                "notes": "",
            }
            summary_rows.append(row)
            print(
                f"  mRNA={len(models)} present={len(present)} frag={len(frag)} weak={len(weak)}"
                + (f" best={best.aa_length}aa id={best.identity:.2f}" if best else "")
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
