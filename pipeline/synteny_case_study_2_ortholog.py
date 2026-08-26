#!/usr/bin/env python3
"""Case study 2: single-gene ortholog synteny across related insect genomes (Inx2, shakB)."""

from __future__ import annotations

import argparse
import csv
import json
import shutil
import sys
from pathlib import Path

from pipeline.common import METADATA_DIR, SYNVOY_DIR, write_csv

CASE_STUDIES_DIR = Path("project/results/synteny_case_studies/case2_ortholog_synteny")
RESULTS_ROOT = SYNVOY_DIR / "results" / "innexin_synvoy"

CASE_GENES = (
    {
        "gene": "Inx2",
        "accession": "Q9V427",
        "job": "dmel_inx2",
        "primary_locus": "synteny_block_locus_1",
        "secondary_locus": "synteny_block_locus_3",
    },
    {
        "gene": "shakB",
        "accession": "P33085",
        "job": "dmel_shakB",
        "primary_locus": "synteny_block_locus_1",
        "secondary_locus": "synteny_block_locus_2",
    },
)


def load_species_map(job: str) -> dict[str, str]:
    path = RESULTS_ROOT / job / "downloaded_genomes/easy_mode_genomes/species_mapping.tsv"
    mapping: dict[str, str] = {}
    if not path.exists():
        return mapping
    for line in path.read_text(encoding="utf-8").splitlines():
        if not line.strip() or line.startswith("#"):
            continue
        parts = line.split("\t")
        if len(parts) >= 2:
            acc = parts[0].replace(".fna", "")
            mapping[acc] = parts[1]
    return mapping


def load_report(job: str) -> dict | None:
    path = RESULTS_ROOT / job / "synvoy_report.json"
    if not path.exists():
        return None
    return json.loads(path.read_text(encoding="utf-8"))


def classify_support(row: dict[str, str]) -> str:
    confident = int(row.get("confident_goi") or 0)
    probable = int(row.get("probable_goi") or 0)
    goi = int(row.get("goi_annotations") or 0)
    flank = int(row.get("flanking_annotations") or 0)
    if confident > 0:
        return "confident_ortholog"
    if probable > 0:
        return "probable_ortholog"
    if goi > 0:
        return "weak_goi_signal"
    if flank > 0:
        return "flanking_synteny_only"
    return "no_support"


def build_ortholog_rows(gene_meta: dict) -> list[dict[str, str]]:
    report = load_report(gene_meta["job"])
    if not report:
        return []
    species_map = load_species_map(gene_meta["job"])
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
            "query_gene": gene_meta["gene"],
            "accession": gene_meta["accession"],
            "target_genome": genome_key,
            "target_species": species_map.get(genome_key, ""),
            "genes_discovered": str(genes_found),
            "total_annotations": str(ann.get("total_annotations", "")),
            "flanking_annotations": str(ann.get("role_counts", {}).get("flanking", "")),
            "goi_annotations": str(ann.get("goi_annotations", "")),
            "confident_goi": str(class_counts.get("confident_goi", 0)),
            "probable_goi": str(class_counts.get("probable_goi", 0)),
            "ambiguous_goi": str(class_counts.get("ambiguous_goi_family_member", 0)),
            "tandem_goi_copies": str(class_counts.get("tandem_goi_copy", 0)),
            "goi_high_confidence": str(ann.get("goi_confidence_counts", {}).get("HIGH", 0)),
            "goi_medium_confidence": str(ann.get("goi_confidence_counts", {}).get("MEDIUM", 0)),
            "goi_low_confidence": str(ann.get("goi_confidence_counts", {}).get("LOW", 0)),
            "synteny_regions": str(reg.get("total_regions", "")),
            "goi_anchor_regions": str(reg.get("goi_anchor_regions", "")),
            "best_region_score": str(reg.get("best_score", "")),
            "ortholog_support": "",
        }
        row["ortholog_support"] = classify_support(row)
        rows.append(row)
    return rows


def build_gene_summary(ortholog_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    by_gene: dict[str, list[dict[str, str]]] = {}
    for row in ortholog_rows:
        by_gene.setdefault(row["query_gene"], []).append(row)

    summaries: list[dict[str, str]] = []
    for gene, rows in sorted(by_gene.items()):
        support_counts: dict[str, int] = {}
        for row in rows:
            support_counts[row["ortholog_support"]] = support_counts.get(row["ortholog_support"], 0) + 1
        best = sorted(
            rows,
            key=lambda r: (
                int(r.get("confident_goi") or 0),
                int(r.get("probable_goi") or 0),
                int(r.get("goi_annotations") or 0),
                float(r.get("best_region_score") or 0),
            ),
            reverse=True,
        )[0]
        summaries.append(
            {
                "query_gene": gene,
                "genomes_searched": str(len(rows)),
                "confident_ortholog_species": str(support_counts.get("confident_ortholog", 0)),
                "probable_ortholog_species": str(support_counts.get("probable_ortholog", 0)),
                "weak_goi_species": str(support_counts.get("weak_goi_signal", 0)),
                "flanking_only_species": str(support_counts.get("flanking_synteny_only", 0)),
                "best_supported_species": best["target_species"],
                "best_supported_genome": best["target_genome"],
                "best_support_class": best["ortholog_support"],
            }
        )
    return summaries


def build_inx2_vs_shakb(ortholog_rows: list[dict[str, str]]) -> list[dict[str, str]]:
    inx2 = {r["target_genome"]: r for r in ortholog_rows if r["query_gene"] == "Inx2"}
    shakb = {r["target_genome"]: r for r in ortholog_rows if r["query_gene"] == "shakB"}
    shared_genomes = sorted(set(inx2) & set(shakb))

    rows: list[dict[str, str]] = []
    for genome in shared_genomes:
        a, b = inx2[genome], shakb[genome]
        rows.append(
            {
                "target_genome": genome,
                "target_species": a["target_species"] or b["target_species"],
                "Inx2_support": a["ortholog_support"],
                "shakB_support": b["ortholog_support"],
                "Inx2_goi_count": a["goi_annotations"],
                "shakB_goi_count": b["goi_annotations"],
                "Inx2_confident_goi": a["confident_goi"],
                "shakB_confident_goi": b["confident_goi"],
                "pattern": f"Inx2={a['ortholog_support']} | shakB={b['ortholog_support']}",
            }
        )
    return rows


def export_figures(out_dir: Path) -> list[dict[str, str]]:
    fig_dir = out_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    rows: list[dict[str, str]] = []
    for meta in CASE_GENES:
        for locus_key in ("primary_locus", "secondary_locus"):
            locus = meta[locus_key]
            src = RESULTS_ROOT / meta["job"] / f"{locus}_synteny_plot.html"
            if not src.exists():
                continue
            dst_name = f"{meta['gene']}_{locus}_synteny_plot.html"
            shutil.copy2(src, fig_dir / dst_name)
            rows.append(
                {
                    "query_gene": meta["gene"],
                    "locus": locus,
                    "exported_figure": str((fig_dir / dst_name).as_posix()),
                    "source_plot": str(src.relative_to(SYNVOY_DIR)),
                }
            )
    return rows


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Case study 2: ortholog synteny across related genomes.")
    parser.add_argument("--output-dir", type=Path, default=CASE_STUDIES_DIR)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    ortholog_rows: list[dict[str, str]] = []
    for meta in CASE_GENES:
        ortholog_rows.extend(build_ortholog_rows(meta))

    if not ortholog_rows:
        print("No SynVoy ortholog data found for Inx2 or shakB.", file=sys.stderr)
        return 1

    summary = build_gene_summary(ortholog_rows)
    comparison = build_inx2_vs_shakb(ortholog_rows)
    figures = export_figures(args.output_dir)

    write_csv(args.output_dir / "case2_ortholog_support.csv", ortholog_rows)
    write_csv(args.output_dir / "case2_gene_summary.csv", summary)
    write_csv(args.output_dir / "case2_inx2_vs_shakB.csv", comparison)
    write_csv(args.output_dir / "case2_figure_index.csv", figures)

    lines = [
        "# Case Study 2: Ortholog synteny across related genomes",
        "",
        "Compares **Inx2** and **shakB** from *Drosophila melanogaster* against SynVoy-selected related genomes.",
        "",
        "## Target species searched",
        "- *Drosophila helvetica*, *D. phalerata*, *D. busckii*, *D. limbata* (Drosophilidae)",
        "- *Anopheles rivulorum* (Diptera; more distant; Inx2 only)",
        "",
        "## Support classes",
        "- `confident_ortholog` — high-confidence GOI",
        "- `probable_ortholog` — probable GOI",
        "- `weak_goi_signal` — GOI candidates, ambiguous only",
        "- `flanking_synteny_only` — conserved neighborhood, no GOI",
        "",
        "## Summary",
    ]
    for row in summary:
        lines.append(
            f"- **{row['query_gene']}**: searched {row['genomes_searched']} genomes; "
            f"best support in *{row['best_supported_species']}* ({row['best_support_class']})"
        )
    lines.extend(
        [
            "",
            "## Tables",
            "- `case2_ortholog_support.csv`",
            "- `case2_gene_summary.csv`",
            "- `case2_inx2_vs_shakB.csv`",
            "",
            "## Figures",
            "Multi-species synteny HTML plots in `figures/`.",
        ]
    )
    (args.output_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")

    print(f"Case study 2 written -> {args.output_dir}")
    for row in summary:
        print(f"  {row['query_gene']}: best={row['best_supported_species']} ({row['best_support_class']})")
    return 0


if __name__ == "__main__":
    sys.exit(main())
