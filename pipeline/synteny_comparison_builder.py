#!/usr/bin/env python3
"""Summarize and compare SynVoy synteny case-study results across query genes."""

from __future__ import annotations

import argparse
import csv
import json
import sys
from pathlib import Path

from pipeline.common import METADATA_DIR, SYNVOY_DIR, write_csv

DEFAULT_QUEUE = METADATA_DIR / "innexin_synvoy_queue.tsv"
DEFAULT_OUT_DIR = METADATA_DIR
RESULTS_ROOT = SYNVOY_DIR / "results" / "innexin_synvoy"


def load_queue(path: Path) -> dict[str, dict[str, str]]:
    lookup: dict[str, dict[str, str]] = {}
    with path.open("r", encoding="utf-8", newline="") as handle:
        for row in csv.DictReader(handle, delimiter="\t"):
            lookup[row["outdir_name"]] = row
    return lookup


def count_synteny_blocks(job_dir: Path) -> int:
    return len(list(job_dir.glob("synteny_block_*_synteny_plot.html")))


def parse_home_flanking(bed_path: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    if not bed_path.exists():
        return rows
    for line in bed_path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) < 4:
            continue
        chrom, start, end, name = parts[0], parts[1], parts[2], parts[3]
        strand = parts[5] if len(parts) > 5 else ""
        is_goi = name.startswith("GOI_") or "GOI" in name.upper()
        rows.append(
            {
                "chrom": chrom,
                "start": start,
                "end": end,
                "gene_id": name,
                "strand": strand,
                "is_goi": "yes" if is_goi else "no",
            }
        )
    return rows


def load_report(job_dir: Path) -> dict | None:
    report_path = job_dir / "synvoy_report.json"
    if not report_path.exists():
        return None
    return json.loads(report_path.read_text(encoding="utf-8"))


def iter_jobs(queue: dict[str, dict[str, str]]) -> list[tuple[dict[str, str], Path]]:
    jobs: list[tuple[dict[str, str], Path]] = []
    if not RESULTS_ROOT.exists():
        return jobs
    for outdir_name, meta in queue.items():
        job_dir = RESULTS_ROOT / outdir_name
        if job_dir.exists():
            jobs.append((meta, job_dir))
    return jobs


def build_job_summary(queue_path: Path) -> list[dict[str, str]]:
    queue = load_queue(queue_path)
    rows: list[dict[str, str]] = []
    for meta, job_dir in iter_jobs(queue):
        report = load_report(job_dir)
        summary = (report or {}).get("summary", {})
        rows.append(
            {
                "gene": meta["gene"],
                "accession": meta["accession"],
                "source_species": meta.get("source_species", ""),
                "outdir_name": meta["outdir_name"],
                "status": "complete" if count_synteny_blocks(job_dir) else "incomplete",
                "synteny_block_plots": str(count_synteny_blocks(job_dir)),
                "genomes_searched": str((report or {}).get("qc_summary", {}).get("total_genomes", "")),
                "genomes_with_hits": str(summary.get("genomes_with_hits", "")),
                "total_hits": str(summary.get("total_hits", "")),
                "genomes_with_annotations": str(summary.get("genomes_with_annotations", "")),
                "total_goi_annotations": str(summary.get("total_goi_annotations", "")),
                "confident_goi_annotations": str(
                    sum(
                        row.get("goi_class_counts", {}).get("confident_goi", 0)
                        for row in (report or {}).get("annotations", {}).get("per_genome", [])
                    )
                ),
                "low_confidence_regions": str(summary.get("low_confidence_regions", "")),
                "report_path": str(job_dir / "synvoy_report.json") if report else "",
            }
        )
    rows.sort(key=lambda r: (r["status"] != "complete", r["gene"].lower()))
    return rows


def build_genome_comparison(queue_path: Path) -> list[dict[str, str]]:
    queue = load_queue(queue_path)
    rows: list[dict[str, str]] = []
    for meta, job_dir in iter_jobs(queue):
        report = load_report(job_dir)
        if not report:
            continue
        genes_discovered = report.get("synteny_results", {}).get("genes_discovered", {})
        hits = report.get("synteny_results", {}).get("synteny_hits_count", {})
        per_genome = {
            row["genome"]: row for row in report.get("annotations", {}).get("per_genome", [])
        }
        genomes = sorted(set(genes_discovered) | set(hits) | set(per_genome))
        for genome in genomes:
            ann = per_genome.get(genome, {})
            rows.append(
                {
                    "gene": meta["gene"],
                    "accession": meta["accession"],
                    "source_species": meta.get("source_species", ""),
                    "target_genome": genome,
                    "genes_discovered": str(genes_discovered.get(genome, "")),
                    "synteny_hits": str(hits.get(genome, "")),
                    "total_annotations": str(ann.get("total_annotations", "")),
                    "flanking_annotations": str(ann.get("role_counts", {}).get("flanking", "")),
                    "goi_annotations": str(ann.get("goi_annotations", "")),
                    "confident_goi": str(ann.get("goi_class_counts", {}).get("confident_goi", "")),
                    "probable_goi": str(ann.get("goi_class_counts", {}).get("probable_goi", "")),
                    "ambiguous_goi": str(ann.get("goi_class_counts", {}).get("ambiguous_goi_family_member", "")),
                    "goi_high_confidence": str(ann.get("goi_confidence_counts", {}).get("HIGH", "")),
                    "goi_medium_confidence": str(ann.get("goi_confidence_counts", {}).get("MEDIUM", "")),
                    "goi_low_confidence": str(ann.get("goi_confidence_counts", {}).get("LOW", "")),
                }
            )
    return rows


def build_home_flanking_table(queue_path: Path) -> list[dict[str, str]]:
    queue = load_queue(queue_path)
    rows: list[dict[str, str]] = []
    for meta, job_dir in iter_jobs(queue):
        for bed_path in sorted(job_dir.glob("plot_inputs_synteny_block_locus_*/synteny_block_locus_*.bed")):
            locus = bed_path.parent.name.replace("plot_inputs_", "")
            genes = parse_home_flanking(bed_path)
            goi_rows = [g for g in genes if g["is_goi"] == "yes"]
            flank_rows = [g for g in genes if g["is_goi"] == "no"]
            flank_names = [g["gene_id"] for g in flank_rows]
            rows.append(
                {
                    "gene": meta["gene"],
                    "accession": meta["accession"],
                    "source_species": meta.get("source_species", ""),
                    "locus": locus,
                    "chrom": genes[0]["chrom"] if genes else "",
                    "flanking_gene_count": str(len(flank_rows)),
                    "goi_count": str(len(goi_rows)),
                    "flanking_genes": ";".join(flank_names[:20]),
                    "goi_ids": ";".join(g["gene_id"] for g in goi_rows),
                    "bed_path": str(bed_path.relative_to(SYNVOY_DIR)),
                }
            )
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Compare SynVoy synteny case-study results.")
    parser.add_argument("--queue", type=Path, default=DEFAULT_QUEUE)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUT_DIR)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if not args.queue.exists():
        print(f"Queue not found: {args.queue}", file=sys.stderr)
        return 1

    job_summary = build_job_summary(args.queue)
    genome_comparison = build_genome_comparison(args.queue)
    home_flanking = build_home_flanking_table(args.queue)

    out_job = args.output_dir / "synteny_job_summary.csv"
    out_genome = args.output_dir / "synteny_genome_comparison.csv"
    out_flank = args.output_dir / "synteny_home_flanking_genes.csv"
    write_csv(out_job, job_summary)
    write_csv(out_genome, genome_comparison)
    write_csv(out_flank, home_flanking)

    complete = sum(1 for row in job_summary if row["status"] == "complete")
    print(f"Wrote {len(job_summary)} job summaries ({complete} complete) -> {out_job}")
    print(f"Wrote {len(genome_comparison)} genome comparisons -> {out_genome}")
    print(f"Wrote {len(home_flanking)} home-locus flanking tables -> {out_flank}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
