#!/usr/bin/env python3
"""Case study 1: Drosophila innexin cluster — ogre vs Inx2 vs Inx3 flanking-gene comparison."""

from __future__ import annotations

import argparse
import csv
import shutil
import sys
from pathlib import Path

from pipeline.common import METADATA_DIR, SYNVOY_DIR, write_csv
from pipeline.synteny_comparison_builder import parse_home_flanking

CASE_STUDIES_DIR = Path("project/results/synteny_case_studies/case1_dmel_innexin_cluster")
RESULTS_ROOT = SYNVOY_DIR / "results" / "innexin_synvoy"

# Comparable SynVoy loci on Dmel X (NC_004354.4)
REGION_A = {
    "label": "proximal_cluster_6.8-7.1Mb",
    "chrom": "NC_004354.4",
    "start": 6_800_000,
    "end": 7_150_000,
    "blocks": [
        ("ogre", "dmel_ogre_inx1", "synteny_block_locus_1"),
        ("Inx2", "dmel_inx2", "synteny_block_locus_1"),
    ],
}
REGION_B = {
    "label": "distal_cluster_20.5-21.0Mb",
    "chrom": "NC_004354.4",
    "start": 20_500_000,
    "end": 21_050_000,
    "blocks": [
        ("ogre", "dmel_ogre_inx1", "synteny_block_locus_2"),
        ("Inx2", "dmel_inx2", "synteny_block_locus_3"),
        ("Inx3", "dmel_inx3", "synteny_block_locus_2"),
    ],
}


def bed_path(job: str, locus: str) -> Path:
    return RESULTS_ROOT / job / f"plot_inputs_{locus}" / f"{locus}.bed"


def plot_path(job: str, locus: str) -> Path:
    return RESULTS_ROOT / job / f"{locus}_synteny_plot.html"


def flanking_gene_set(bed: Path) -> set[str]:
    return {row["gene_id"] for row in parse_home_flanking(bed) if row["is_goi"] == "no"}


def goi_gene_set(bed: Path) -> set[str]:
    return {row["gene_id"] for row in parse_home_flanking(bed) if row["is_goi"] == "yes"}


def compare_region(region: dict, all_queries: list[str]) -> tuple[list[dict[str, str]], list[dict[str, str]]]:
    genes_by_query: dict[str, set[str]] = {}
    goi_by_query: dict[str, set[str]] = {}
    for gene, job, locus in region["blocks"]:
        path = bed_path(job, locus)
        if not path.exists():
            continue
        genes_by_query[gene] = flanking_gene_set(path)
        goi_by_query[gene] = goi_gene_set(path)

    if not genes_by_query:
        return [], []

    all_genes = sorted(set().union(*genes_by_query.values()))
    queries = list(genes_by_query.keys())
    matrix_rows: list[dict[str, str]] = []
    for gene_id in all_genes:
        row = {
            "region": region["label"],
            "flanking_gene_id": gene_id,
        }
        present_in = []
        for query in all_queries:
            present = gene_id in genes_by_query.get(query, set())
            row[query] = "yes" if present else "no"
            if present:
                present_in.append(query)
        row["present_in_count"] = str(len(present_in))
        row["present_in_genes"] = ";".join(present_in)
        if len(present_in) == len(queries):
            row["sharing_category"] = "shared_all"
        elif len(present_in) == 1:
            row["sharing_category"] = f"unique_{present_in[0]}"
        elif len(present_in) == 0:
            row["sharing_category"] = "none"
        else:
            row["sharing_category"] = "partial_shared"
        matrix_rows.append(row)

    summary_rows: list[dict[str, str]] = []
    shared_all = {g for g in all_genes if all(g in genes_by_query[q] for q in queries)}
    for query in queries:
        unique = genes_by_query[query] - set().union(*(genes_by_query[q] for q in queries if q != query))
        summary_rows.append(
            {
                "region": region["label"],
                "query_gene": query,
                "flanking_gene_count": str(len(genes_by_query[query])),
                "goi_model_count": str(len(goi_by_query.get(query, set()))),
                "shared_with_all_paralogs": str(len(shared_all)),
                "unique_to_query": str(len(unique)),
                "shared_gene_ids": ";".join(sorted(shared_all)[:30]) + ("..." if len(shared_all) > 30 else ""),
                "unique_gene_ids": ";".join(sorted(unique)[:30]) + ("..." if len(unique) > 30 else ""),
            }
        )
    return matrix_rows, summary_rows


def export_figure_links(out_dir: Path) -> list[dict[str, str]]:
    rows: list[dict[str, str]] = []
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    for region, blocks in ((REGION_A, REGION_A["blocks"]), (REGION_B, REGION_B["blocks"])):
        for gene, job, locus in blocks:
            src = plot_path(job, locus)
            if not src.exists():
                continue
            dst_name = f"{gene}_{locus}_synteny_plot.html"
            shutil.copy2(src, fig_dir / dst_name)
            rows.append(
                {
                    "region": region["label"],
                    "query_gene": gene,
                    "locus": locus,
                    "source_plot": str(src.relative_to(SYNVOY_DIR)),
                    "exported_figure": str((fig_dir / dst_name).as_posix()),
                }
            )
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Case study 1: Dmel innexin cluster synteny comparison.")
    parser.add_argument("--output-dir", type=Path, default=CASE_STUDIES_DIR)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    matrix: list[dict[str, str]] = []
    summary: list[dict[str, str]] = []
    all_queries = sorted({gene for region in (REGION_A, REGION_B) for gene, _, _ in region["blocks"]})
    for region in (REGION_A, REGION_B):
        region_matrix, region_summary = compare_region(region, all_queries)
        matrix.extend(region_matrix)
        summary.extend(region_summary)

    if not matrix:
        print("No comparable SynVoy blocks found. Run ogre, Inx2, Inx3 first.", file=sys.stderr)
        return 1

    figure_links = export_figure_links(args.output_dir)

    write_csv(args.output_dir / "case1_flanking_gene_matrix.csv", matrix)
    write_csv(args.output_dir / "case1_region_summary.csv", summary)
    write_csv(args.output_dir / "case1_figure_index.csv", figure_links)

    # Short narrative for thesis notes
    lines = [
        "# Case Study 1: Drosophila innexin cluster synteny",
        "",
        "Compares home-genome flanking genes around **ogre**, **Inx2**, and **Inx3**.",
        "",
        "## Regions",
        f"- **Proximal cluster** ({REGION_A['label']}): ogre locus 1 vs Inx2 locus 1",
        f"- **Distal cluster** ({REGION_B['label']}): ogre locus 2 vs Inx2 locus 3 vs Inx3 locus 2",
        "",
        "## Summary",
    ]
    for row in summary:
        lines.append(
            f"- **{row['query_gene']}** @ {row['region']}: "
            f"{row['flanking_gene_count']} flanking genes, "
            f"{row['shared_with_all_paralogs']} shared with all compared paralogs, "
            f"{row['unique_to_query']} unique"
        )
    lines.extend(
        [
            "",
            "## Figures",
            "Open HTML files in `figures/` or use `case1_figure_index.csv`.",
            "",
            "## Tables",
            "- `case1_flanking_gene_matrix.csv` — presence/absence per flanking gene",
            "- `case1_region_summary.csv` — shared vs unique counts",
        ]
    )
    (args.output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Case study 1 written -> {args.output_dir}")
    for row in summary:
        print(
            f"  {row['query_gene']} @ {row['region']}: "
            f"shared={row['shared_with_all_paralogs']} unique={row['unique_to_query']}"
        )
    return 0


if __name__ == "__main__":
    sys.exit(main())
