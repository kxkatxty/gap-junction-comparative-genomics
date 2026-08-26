#!/usr/bin/env python3
"""Case study 3: cross-clade synteny — C. elegans inx-2 vs Drosophila innexins."""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from pathlib import Path

from pipeline.common import SYNVOY_DIR, write_csv
from pipeline.synteny_case_study_2_ortholog import (
    CASE_GENES as FLY_CASE_GENES,
    build_ortholog_rows,
    build_gene_summary,
    classify_support,
    load_report,
    load_species_map,
)
from pipeline.synteny_comparison_builder import parse_home_flanking

CASE_STUDIES_DIR = Path("project/results/synteny_case_studies/case3_cross_clade_inx2")
RESULTS_ROOT = SYNVOY_DIR / "results" / "innexin_synvoy"

CELE_INX2 = {
    "gene": "inx-2",
    "accession": "Q9U3K5",
    "job": "cele_inx-2",
    "clade": "nematode",
    "source_species": "Caenorhabditis elegans",
}
FLY_REFERENCE = {
    "gene": "Inx2",
    "accession": "Q9V427",
    "job": "dmel_inx2",
    "clade": "insect",
    "source_species": "Drosophila melanogaster",
}


def bed_path(job: str, locus: str = "synteny_block_locus_1") -> Path:
    return RESULTS_ROOT / job / f"plot_inputs_{locus}" / f"{locus}.bed"


def home_locus_stats(job: str) -> dict[str, str]:
    path = bed_path(job)
    if not path.exists():
        return {}
    genes = parse_home_flanking(path)
    flank = [g for g in genes if g["is_goi"] == "no"]
    goi = [g for g in genes if g["is_goi"] == "yes"]
    return {
        "home_flanking_count": str(len(flank)),
        "home_goi_count": str(len(goi)),
        "home_chrom": genes[0]["chrom"] if genes else "",
    }


def build_cele_ortholog_rows() -> list[dict[str, str]]:
    report = load_report(CELE_INX2["job"])
    if not report:
        return []
    species_map = load_species_map(CELE_INX2["job"])
    per_genome = {row["genome"].replace(".fna", ""): row for row in report.get("annotations", {}).get("per_genome", [])}
    regions = {row["genome"].replace(".fna", ""): row for row in report.get("regions", {}).get("per_genome", [])}
    discovered = report.get("synteny_results", {}).get("genes_discovered", {})

    rows: list[dict[str, str]] = []
    for genome_acc, genes_found in sorted(discovered.items()):
        genome_key = genome_acc.replace(".fna", "")
        ann = per_genome.get(genome_key, {})
        reg = regions.get(genome_key, {})
        class_counts = ann.get("goi_class_counts", {})
        row = {
            "query_gene": CELE_INX2["gene"],
            "query_clade": CELE_INX2["clade"],
            "accession": CELE_INX2["accession"],
            "source_species": CELE_INX2["source_species"],
            "target_genome": genome_key,
            "target_species": species_map.get(genome_key, ""),
            "genes_discovered": str(genes_found),
            "goi_annotations": str(ann.get("goi_annotations", "")),
            "confident_goi": str(class_counts.get("confident_goi", 0)),
            "probable_goi": str(class_counts.get("probable_goi", 0)),
            "ambiguous_goi": str(class_counts.get("ambiguous_goi_family_member", 0)),
            "flanking_annotations": str(ann.get("role_counts", {}).get("flanking", "")),
            "synteny_regions": str(reg.get("total_regions", "")),
            "best_region_score": str(reg.get("best_score", "")),
            "ortholog_support": "",
        }
        row["ortholog_support"] = classify_support(row)
        rows.append(row)
    return rows


def build_clade_comparison() -> list[dict[str, str]]:
    cele_rows = build_cele_ortholog_rows()
    fly_rows = [r for r in build_ortholog_rows(FLY_REFERENCE) if r["query_gene"] == "Inx2"]
    if not cele_rows and not fly_rows:
        return []

    def summarize(clade: str, gene: str, source: str, rows: list[dict[str, str]], home: dict[str, str]) -> dict[str, str]:
        support: dict[str, int] = {}
        for row in rows:
            support[row["ortholog_support"]] = support.get(row["ortholog_support"], 0) + 1
        return {
            "clade": clade,
            "query_gene": gene,
            "source_species": source,
            "genomes_searched": str(len(rows)),
            "confident_ortholog_species": str(support.get("confident_ortholog", 0)),
            "probable_ortholog_species": str(support.get("probable_ortholog", 0)),
            "weak_goi_species": str(support.get("weak_goi_signal", 0)),
            "flanking_only_species": str(support.get("flanking_synteny_only", 0)),
            "total_goi_annotations": str(sum(int(r.get("goi_annotations") or 0) for r in rows)),
            "home_flanking_count": home.get("home_flanking_count", ""),
            "home_goi_count": home.get("home_goi_count", ""),
        }

    return [
        summarize(
            CELE_INX2["clade"],
            CELE_INX2["gene"],
            CELE_INX2["source_species"],
            cele_rows,
            home_locus_stats(CELE_INX2["job"]),
        ),
        summarize(
            FLY_REFERENCE["clade"],
            FLY_REFERENCE["gene"],
            FLY_REFERENCE["source_species"],
            fly_rows,
            home_locus_stats(FLY_REFERENCE["job"]),
        ),
    ]


def export_figures(out_dir: Path) -> list[dict[str, str]]:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    for locus in sorted((RESULTS_ROOT / CELE_INX2["job"]).glob("synteny_block_*_synteny_plot.html")):
        dst = fig_dir / f"cele_inx-2_{locus.name}"
        shutil.copy2(locus, dst)
        rows.append({"query_gene": "inx-2", "figure": str(dst.as_posix())})
    fly_plot = RESULTS_ROOT / FLY_REFERENCE["job"] / "synteny_block_locus_1_synteny_plot.html"
    if fly_plot.exists():
        dst = fig_dir / "dmel_Inx2_synteny_block_locus_1_synteny_plot.html"
        shutil.copy2(fly_plot, dst)
        rows.append({"query_gene": "Inx2", "figure": str(dst.as_posix())})
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Case study 3: cross-clade inx-2 vs fly Inx2 synteny.")
    parser.add_argument("--output-dir", type=Path, default=CASE_STUDIES_DIR)
    parser.add_argument("--allow-incomplete", action="store_true", help="Write fly-only comparison if cele run pending.")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    cele_rows = build_cele_ortholog_rows()
    if not cele_rows and not args.allow_incomplete:
        print(
            "C. elegans inx-2 SynVoy results not found. Run:\n"
            "  python3 -m pipeline.case_study_synteny_runner --gene inx-2",
            file=sys.stderr,
        )
        return 1

    comparison = build_clade_comparison()
    figures = export_figures(args.output_dir)

    write_csv(args.output_dir / "case3_cele_inx2_ortholog_support.csv", cele_rows)
    write_csv(args.output_dir / "case3_cross_clade_comparison.csv", comparison)
    write_csv(args.output_dir / "case3_figure_index.csv", figures)

    lines = [
        "# Case Study 3: Cross-clade synteny (C. elegans inx-2 vs Drosophila Inx2)",
        "",
        "Compares ortholog/synteny support for the same innexin family across **nematode** vs **insect** queries.",
        "",
    ]
    for row in comparison:
        lines.append(
            f"- **{row['source_species']} {row['query_gene']}** ({row['clade']}): "
            f"{row['genomes_searched']} genomes searched; "
            f"confident={row['confident_ortholog_species']}, "
            f"probable={row['probable_ortholog_species']}, "
            f"weak={row['weak_goi_species']}, "
            f"flanking-only={row['flanking_only_species']}"
        )
    if not cele_rows:
        lines.append("\n> C. elegans inx-2 SynVoy run pending — rerun after `case_study_synteny_runner --gene inx-2`.")
    (args.output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Case study 3 written -> {args.output_dir}")
    for row in comparison:
        print(f"  {row['clade']} {row['query_gene']}: genomes={row['genomes_searched']} GOI total={row['total_goi_annotations']}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
