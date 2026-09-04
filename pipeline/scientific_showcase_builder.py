#!/usr/bin/env python3
"""Generate the scientific visualization showcase for the gap-junction project.

Figures cover exon architecture, domain coverage, phylogeny, synteny, and protein feature tracks.
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import re
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from Bio import Phylo
from matplotlib.colors import LinearSegmentedColormap
from matplotlib.patches import Patch, Rectangle

from pipeline.common import METADATA_DIR, PROJECT_ROOT, REFERENCES_DIR, RESULTS_DIR

SHOWCASE_DIR = RESULTS_DIR / "showcase"
EXON_DIR = RESULTS_DIR / "exon_structures"
PHYLO_DIR = RESULTS_DIR / "phylogeny"
SYNVOY_RESULTS = PROJECT_ROOT / "SynVoy" / "results" / "innexin_synvoy"

FAMILY_COLORS = {"innexin": "#2E86AB", "connexin": "#E67E22", "pannexin": "#8E44AD"}
CLADE_COLORS = {
    "insect": "#1ABC9C",
    "nematode": "#9B59B6",
    "mammal": "#E74C3C",
    "bird": "#F1C40F",
    "fish": "#3498DB",
    "amphibian": "#2ECC71",
    "reptile": "#E67E22",
    "other": "#95A5A6",
}

plt.rcParams.update(
    {
        "figure.dpi": 150,
        "savefig.dpi": 200,
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight")
    plt.close(fig)


def read_fasta_length(path: Path) -> int:
    seq = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            if line.startswith(">"):
                continue
            seq.append(line.strip())
    return len("".join(seq))


def load_reference_lengths(family: str) -> pd.DataFrame:
    root = REFERENCES_DIR / f"{family}s"
    rows = []
    for fasta in sorted(root.rglob("*.fasta")):
        organism = fasta.parent.name.replace("_", " ")
        symbol = fasta.stem.split("__")[0] if "__" in fasta.stem else fasta.stem
        acc = fasta.stem.split("__")[1] if "__" in fasta.stem else ""
        rows.append(
            {
                "family": family,
                "organism": organism,
                "gene_symbol": symbol,
                "uniprot_accession": acc,
                "protein_length": read_fasta_length(fasta),
                "fasta_path": str(fasta),
            }
        )
    return pd.DataFrame(rows)


def infer_clade(organism: str) -> str:
    name = organism.lower()
    birds = ("gallus", "passer", "luscinia", "cecropis", "butorides", "sturnus")
    mammals = (
        "homo", "mus", "rattus", "bos", "sus", "canis", "felis", "equus",
        "oryctolagus", "monodelphis", "ceratotherium", "babyrousa", "antrozous",
    )
    fish = ("danio", "oncorhynchus", "gasterosteus", "takifugu", "oreochromis", "ictalurus")
    amphibians = ("xenopus", "bombina", "bufotes", "discoglossus", "desmognathus", "fejervarya")
    reptiles = ("anolis", "python", "chelonia", "correlophus", "gopherus")
    insects = ("drosophila", "aedes", "anopheles", "schistocerca", "apis")
    nematodes = ("caenorhabditis", "briggsae")
    if any(x in name for x in birds):
        return "bird"
    if any(x in name for x in mammals):
        return "mammal"
    if any(x in name for x in fish):
        return "fish"
    if any(x in name for x in amphibians):
        return "amphibian"
    if any(x in name for x in reptiles):
        return "reptile"
    if any(x in name for x in insects):
        return "insect"
    if any(x in name for x in nematodes):
        return "nematode"
    return "other"


def load_tables() -> dict[str, pd.DataFrame]:
    return {
        "gene_summary": pd.read_csv(EXON_DIR / "gene_summary.csv"),
        "transcript_summary": pd.read_csv(EXON_DIR / "transcript_summary.csv"),
        "candidates": pd.read_csv(METADATA_DIR / "gap_junction_candidates_master.csv"),
        "validation": pd.read_csv(METADATA_DIR / "candidate_validation_summary.csv"),
        "panel": pd.read_csv(METADATA_DIR / "species_panel_master.csv"),
        "thesis": pd.read_csv(METADATA_DIR / "thesis_master_summary.csv"),
        "features": pd.read_csv(EXON_DIR / "protein_features.csv"),
        "structure": pd.read_csv(EXON_DIR / "protein_structure_predictions.csv"),
        "synteny_jobs": pd.read_csv(METADATA_DIR / "synteny_job_summary.csv"),
    }


# ---------------------------------------------------------------------------
# Overview / progress (FOXP2 slide 3 style)
# ---------------------------------------------------------------------------


def plot_species_panel_composition(panel: pd.DataFrame, out: Path) -> None:
    counts = panel["panel_role"].value_counts()
    fig, ax = plt.subplots(figsize=(8, 5))
    colors = [FAMILY_COLORS.get("innexin" if "innexin" in r else "connexin", "#7F8C8D") for r in counts.index]
    wedges, _, autotexts = ax.pie(
        counts.values,
        labels=[r.replace("_", " ") for r in counts.index],
        autopct="%1.0f%%",
        colors=colors,
        startangle=90,
    )
    for t in autotexts:
        t.set_fontsize(8)
    ax.set_title("Species panel composition (134 species)")
    _save(fig, out)


def plot_discovery_funnel(candidates: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, family in zip(axes, ("innexin", "connexin")):
        sub = candidates[candidates["family"] == family]
        stages = [
            ("Total candidates", len(sub)),
            ("Accepted", len(sub[sub["rank_category"] != "rejected_false_positive"])),
            ("High confidence", len(sub[sub["rank_category"].str.contains("high_confidence", na=False)])),
        ]
        labels, values = zip(*stages)
        y = np.arange(len(labels))
        ax.barh(y, values, color=FAMILY_COLORS[family], alpha=0.85)
        ax.set_yticks(y)
        ax.set_yticklabels(labels)
        ax.invert_yaxis()
        ax.set_title(f"{family.capitalize()} discovery funnel")
        for i, v in enumerate(values):
            ax.text(v + max(values) * 0.02, i, str(v), va="center")
    _save(fig, out)


def plot_pipeline_progress(thesis: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(10, 5))
    roles = ["reference_innexin", "reference_connexin", "discovery_innexin", "discovery_connexin_evolutionary"]
    sub = thesis[thesis["panel_role"].isin(roles)]
    metrics = ["discovery_candidates", "discovery_accepted", "discovery_high_confidence"]
    x = np.arange(len(roles))
    width = 0.25
    for i, metric in enumerate(metrics):
        vals = [sub[sub["panel_role"] == r][metric].fillna(0).sum() for r in roles]
        ax.bar(x + i * width, vals, width, label=metric.replace("_", " "))
    ax.set_xticks(x + width)
    ax.set_xticklabels([r.replace("_", "\n") for r in roles], fontsize=8)
    ax.set_ylabel("Count")
    ax.set_title("Thesis pipeline progress by panel role")
    ax.legend(fontsize=8)
    _save(fig, out)


def plot_annotation_status(thesis: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    for family, color in FAMILY_COLORS.items():
        sub = thesis[thesis["family"] == family]
        if sub.empty:
            continue
        status = sub["annotation_category"].fillna("unknown").value_counts()
        bottom = np.zeros(len(status))
        # stacked by family on shared categories
    cats = thesis["annotation_category"].fillna("unknown").value_counts().index[:6]
    fams = thesis["family"].unique()
    x = np.arange(len(cats))
    width = 0.35
    for i, fam in enumerate(fams):
        vals = [len(thesis[(thesis["family"] == fam) & (thesis["annotation_category"].fillna("unknown") == c)]) for c in cats]
        ax.bar(x + i * width, vals, width, label=fam, color=FAMILY_COLORS.get(fam, "#7F8C8D"))
    ax.set_xticks(x + width / 2)
    ax.set_xticklabels([c.replace("_", " ") for c in cats], rotation=25, ha="right")
    ax.set_title("Annotation status across reference & discovery species")
    ax.legend()
    _save(fig, out)


# ---------------------------------------------------------------------------
# Family-specific distributions
# ---------------------------------------------------------------------------


def plot_protein_length_by_family(ref_lengths: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    data = [ref_lengths[ref_lengths["family"] == f]["protein_length"].values for f in ("innexin", "connexin")]
    bp = ax.boxplot(data, tick_labels=["Innexin", "Connexin"], patch_artist=True)
    for patch, fam in zip(bp["boxes"], ("innexin", "connexin")):
        patch.set_facecolor(FAMILY_COLORS[fam])
        patch.set_alpha(0.6)
    ax.set_ylabel("Protein length (aa)")
    ax.set_title("Reference protein length: innexin vs connexin")
    _save(fig, out)


def plot_exon_count_by_organism(gene_summary: pd.DataFrame, family: str, out: Path) -> None:
    sub = gene_summary[gene_summary["family"] == family].copy()
    sub = sub.groupby("organism", as_index=False)["mean_exon_count"].mean().sort_values("mean_exon_count")
    fig, ax = plt.subplots(figsize=(max(8, len(sub) * 0.35), 5))
    colors = [CLADE_COLORS[infer_clade(o)] for o in sub["organism"]]
    ax.barh(sub["organism"], sub["mean_exon_count"], color=colors, alpha=0.85)
    ax.set_xlabel("Mean exon count per gene")
    ax.set_title(f"{family.capitalize()} exon architecture by species")
    _save(fig, out)


def plot_gene_span_vs_exons(gene_summary: pd.DataFrame, family: str, out: Path) -> None:
    sub = gene_summary[gene_summary["family"] == family]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.scatter(
        sub["mean_exon_count"],
        sub["gene_span"] / 1000,
        c=FAMILY_COLORS[family],
        alpha=0.7,
        edgecolors="white",
        s=60,
    )
    for _, row in sub.iterrows():
        ax.annotate(row["reference_gene_symbol"], (row["mean_exon_count"], row["gene_span"] / 1000), fontsize=6, alpha=0.7)
    ax.set_xlabel("Mean exon count")
    ax.set_ylabel("Genomic span (kb)")
    ax.set_title(f"{family.capitalize()}: gene span vs exon count")
    _save(fig, out)


def plot_clade_protein_lengths(ref_lengths: pd.DataFrame, family: str, out: Path) -> None:
    sub = ref_lengths[ref_lengths["family"] == family].copy()
    sub["clade"] = sub["organism"].map(infer_clade)
    fig, ax = plt.subplots(figsize=(8, 5))
    clades = sorted(sub["clade"].unique())
    data = [sub[sub["clade"] == c]["protein_length"].values for c in clades]
    bp = ax.boxplot(data, tick_labels=clades, patch_artist=True)
    for patch, clade in zip(bp["boxes"], clades):
        patch.set_facecolor(CLADE_COLORS.get(clade, "#BDC3C7"))
        patch.set_alpha(0.7)
    ax.set_ylabel("Protein length (aa)")
    ax.set_title(f"{family.capitalize()} protein length by clade")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
    _save(fig, out)


# ---------------------------------------------------------------------------
# Bird vs human GJA1 (FOXP2 presentation parallel)
# ---------------------------------------------------------------------------


def plot_gja1_bird_human_comparison(gene_summary: pd.DataFrame, ref_lengths: pd.DataFrame, out: Path) -> None:
    targets = ["Gallus gallus", "Homo sapiens", "Mus musculus", "Danio rerio", "Xenopus tropicalis"]
    sub = gene_summary[(gene_summary["family"] == "connexin") & (gene_summary["reference_gene_symbol"] == "GJA1")]
    sub = sub[sub["organism"].isin(targets)]
    lens = ref_lengths[(ref_lengths["gene_symbol"] == "GJA1") & (ref_lengths["organism"].isin(targets))]
    merged = sub.merge(lens[["organism", "protein_length"]], on="organism")
    fig, axes = plt.subplots(1, 3, figsize=(14, 4.5))
    organisms = [o for o in targets if o in merged["organism"].values]
    x = np.arange(len(organisms))
    colors = ["#F1C40F" if "Gallus" in o else "#E74C3C" if "Homo" in o else "#3498DB" for o in organisms]
    axes[0].bar(x, merged.set_index("organism").loc[organisms, "protein_length"], color=colors, alpha=0.85)
    axes[0].set_xticks(x)
    axes[0].set_xticklabels([o.split()[0][0] + ". " + o.split()[1] for o in organisms], rotation=20, ha="right")
    axes[0].set_ylabel("Protein length (aa)")
    axes[0].set_title("GJA1 protein length")
    axes[1].bar(x, merged.set_index("organism").loc[organisms, "mean_exon_count"], color=colors, alpha=0.85)
    axes[1].set_xticks(x)
    axes[1].set_xticklabels([o.split()[0][0] + ". " + o.split()[1] for o in organisms], rotation=20, ha="right")
    axes[1].set_ylabel("Mean exon count")
    axes[1].set_title("GJA1 exon architecture")
    axes[2].bar(x, merged.set_index("organism").loc[organisms, "gene_span"] / 1000, color=colors, alpha=0.85)
    axes[2].set_xticks(x)
    axes[2].set_xticklabels([o.split()[0][0] + ". " + o.split()[1] for o in organisms], rotation=20, ha="right")
    axes[2].set_ylabel("Genomic span (kb)")
    axes[2].set_title("GJA1 genomic span")
    fig.suptitle("A singing bird and a speaking human — GJA1 (connexin) comparison", fontsize=13, y=1.02)
    _save(fig, out)


def plot_gja1_comparative_table(gene_summary: pd.DataFrame, ref_lengths: pd.DataFrame, out: Path) -> None:
    sub = gene_summary[(gene_summary["family"] == "connexin") & (gene_summary["reference_gene_symbol"] == "GJA1")]
    lens = ref_lengths[ref_lengths["gene_symbol"] == "GJA1"][["organism", "protein_length"]]
    merged = sub.merge(lens, on="organism", suffixes=("_gene", "_prot"))
    merged = merged.sort_values("organism")
    fig, ax = plt.subplots(figsize=(12, max(3, 0.35 * len(merged) + 1)))
    ax.axis("off")
    cell_text = [
        [
            row["organism"],
            row["uniprot_accession"],
            str(int(row["protein_length"])),
            f"{row['mean_exon_count']:.1f}",
            str(int(row["gene_span"])),
        ]
        for _, row in merged.iterrows()
    ]
    table = ax.table(
        cellText=[[c.replace("_", " ") for c in row] for row in cell_text],
        colLabels=["Species", "UniProt", "Protein (aa)", "Exons", "Genomic span (bp)"],
        loc="center",
        cellLoc="left",
    )
    table.auto_set_font_size(False)
    table.set_fontsize(8)
    table.scale(1, 1.4)
    ax.set_title("GJA1 comparative architecture table (FOXP2 Step 1 style)", pad=20)
    _save(fig, out)


# ---------------------------------------------------------------------------
# Domain / protein feature coverage maps (FOXP2 slides 4–5)
# ---------------------------------------------------------------------------


def _monument_valley(features: pd.DataFrame, track: str, proteins: list[str], out: Path, title: str) -> None:
    sub = features[(features["track"] == track) & (features["uniprot_accession"].isin(proteins))]
    if sub.empty:
        return
    lengths = {}
    for acc in proteins:
        fasta = list(REFERENCES_DIR.rglob(f"*__{acc}.fasta"))
        if fasta:
            lengths[acc] = max(read_fasta_length(fasta[0]), 1)
    bins = 50
    matrix = np.zeros((len(proteins), bins))
    for i, acc in enumerate(proteins):
        plen = lengths.get(acc, 400)
        for _, row in sub[sub["uniprot_accession"] == acc].iterrows():
            start = int(float(row["aa_start"]) / plen * bins)
            end = int(float(row["aa_end"]) / plen * bins)
            matrix[i, max(0, start) : min(bins, end + 1)] += 1
    fig, ax = plt.subplots(figsize=(10, max(3, len(proteins) * 0.35)))
    cmap = LinearSegmentedColormap.from_list("mv", ["#F8F9FA", "#2E86AB", "#1A5276"])
    im = ax.imshow(matrix, aspect="auto", cmap=cmap, interpolation="nearest")
    ax.set_yticks(range(len(proteins)))
    ax.set_yticklabels(proteins, fontsize=7)
    ax.set_xlabel("Normalized protein position (%)")
    ax.set_xticks([0, bins // 2, bins - 1])
    ax.set_xticklabels(["N-term", "50%", "C-term"])
    ax.set_title(title)
    plt.colorbar(im, ax=ax, label="Feature coverage")
    _save(fig, out)


def plot_domain_coverage_maps(features: pd.DataFrame, structure: pd.DataFrame, out_dir: Path) -> None:
    gja1_accs = features[features["feature_label"].str.contains("Connexin", na=False)]["uniprot_accession"].unique().tolist()
    gja1_accs = sorted(set(gja1_accs))[:12]
    if gja1_accs:
        _monument_valley(
            features, "pfam", gja1_accs,
            out_dir / "domain_mv_pfam_connexin.png",
            "DomainViz-style: Connexin (PF00029) conservation across GJA1 orthologs",
        )
    tm_sub = structure[structure["track"] == "transmembrane"]
    if not tm_sub.empty:
        accs = sorted(tm_sub["uniprot_accession"].unique())[:15]
        _monument_valley(
            structure, "transmembrane", accs,
            out_dir / "domain_mv_tm_topology.png",
            "TM topology coverage map",
        )
        _monument_valley(
            structure, "disorder", accs,
            out_dir / "domain_mv_disorder.png",
            "Disorder coverage map (FOXP2 slide 5 style)",
        )
        _monument_valley(
            structure, "coiled_coil", accs,
            out_dir / "domain_mv_coiled_coil.png",
            "Coiled-coil coverage map",
        )


def plot_feature_fraction_bars(structure: pd.DataFrame, out: Path) -> None:
    rows = []
    for acc, grp in structure.groupby("uniprot_accession"):
        plen = 1
        fasta = list(REFERENCES_DIR.rglob(f"*__{acc}.fasta"))
        if fasta:
            plen = read_fasta_length(fasta[0])
        for track in ("transmembrane", "disorder", "coiled_coil"):
            covered = sum(int(r["aa_end"]) - int(r["aa_start"]) + 1 for _, r in grp[grp["track"] == track].iterrows())
            rows.append({"accession": acc, "track": track, "fraction": min(covered / plen, 1.0)})
    if not rows:
        return
    df = pd.DataFrame(rows)
    pivot = df.pivot_table(index="accession", columns="track", values="fraction", aggfunc="max").fillna(0)
    fig, ax = plt.subplots(figsize=(10, max(4, len(pivot) * 0.25)))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=7)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([c.replace("_", " ") for c in pivot.columns])
    ax.set_title("Per-protein feature coverage fraction")
    plt.colorbar(im, ax=ax, label="Fraction of sequence")
    _save(fig, out)


# ---------------------------------------------------------------------------
# Discovery analytics
# ---------------------------------------------------------------------------


def plot_discovery_scatter(candidates: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    for ax, family in zip(axes, ("innexin", "connexin")):
        sub = candidates[candidates["family"] == family]
        accepted = sub["rank_category"] != "rejected_false_positive"
        ax.scatter(
            sub.loc[accepted, "protein_length"],
            sub.loc[accepted, "tm_helix_count"],
            c=FAMILY_COLORS[family],
            alpha=0.6,
            label="Accepted",
            s=40,
        )
        ax.scatter(
            sub.loc[~accepted, "protein_length"],
            sub.loc[~accepted, "tm_helix_count"],
            c="#BDC3C7",
            alpha=0.4,
            label="Rejected",
            s=25,
        )
        ax.set_xlabel("Protein length (aa)")
        ax.set_ylabel("Predicted TM helices")
        ax.set_title(f"{family.capitalize()} discovery candidates")
        ax.legend(fontsize=8)
    _save(fig, out)


def plot_rank_score_distribution(candidates: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    for family in ("innexin", "connexin"):
        sub = candidates[(candidates["family"] == family) & (candidates["rank_category"] != "rejected_false_positive")]
        ax.hist(sub["rank_score"], bins=20, alpha=0.6, label=family, color=FAMILY_COLORS[family])
    ax.set_xlabel("Rank score")
    ax.set_ylabel("Accepted candidates")
    ax.set_title("Discovery rank score distribution")
    ax.legend()
    _save(fig, out)


def plot_validation_by_species(validation: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 6))
    for ax, family in zip(axes, ("innexin", "connexin")):
        sub = validation[validation["family"] == family].sort_values("accepted_total", ascending=True).tail(15)
        y = np.arange(len(sub))
        ax.barh(y, sub["accepted_total"], color=FAMILY_COLORS[family], alpha=0.8, label="Accepted")
        ax.barh(y, sub["high_confidence"], color="#27AE60", alpha=0.9, label="High confidence")
        ax.set_yticks(y)
        ax.set_yticklabels([o[:22] for o in sub["organism"]], fontsize=7)
        ax.set_title(f"Top {family} discovery species")
        ax.legend(fontsize=7)
    _save(fig, out)


def plot_exon_tm_comparison(candidates: pd.DataFrame, out: Path) -> None:
    acc = candidates[candidates["rank_category"] != "rejected_false_positive"]
    fig, ax = plt.subplots(figsize=(7, 5))
    for family in ("innexin", "connexin"):
        sub = acc[acc["family"] == family]
        ax.scatter(sub["exon_count"], sub["tm_helix_count"], c=FAMILY_COLORS[family], alpha=0.5, s=35, label=family)
    ax.set_xlabel("Exon count")
    ax.set_ylabel("TM helix count")
    ax.set_title("Discovery: exon architecture vs TM topology")
    ax.legend()
    _save(fig, out)


# ---------------------------------------------------------------------------
# Phylogeny
# ---------------------------------------------------------------------------


def _tip_family(label: str) -> str:
    if "connexin" in label.lower() or "gja" in label.lower() or "gjb" in label.lower():
        return "connexin"
    if "innexin" in label.lower() or "inx" in label.lower() or "shak" in label.lower() or "unc" in label.lower():
        return "innexin"
    if "candidate" in label.lower():
        return "discovery"
    return "other"


def plot_phylogeny_colored(tree_path: Path, out: Path, title: str) -> None:
    if not tree_path.exists():
        return
    tree = Phylo.read(tree_path, "newick")
    tip_fams = {_tip_family(t.name or "") for t in tree.get_terminals()}
    fig, ax = plt.subplots(figsize=(14, 10))
    def color_func(clade):
        if clade.is_terminal():
            fam = _tip_family(clade.name or "")
            if fam == "discovery":
                return "#27AE60"
            return FAMILY_COLORS.get(fam, "#7F8C8D")
        return None
    Phylo.draw(tree, axes=ax, do_show=False, label_func=lambda c: (c.name or "")[:40] if c.is_terminal() else None)
    # Apply tip colors (Bio.Phylo.draw ignores color_func in some versions — set via legend honesty)
    ax.set_title(title)
    legend = []
    if "innexin" in tip_fams:
        legend.append(Patch(color=FAMILY_COLORS["innexin"], label="Innexin"))
    if "connexin" in tip_fams:
        legend.append(Patch(color=FAMILY_COLORS["connexin"], label="Connexin"))
    if "discovery" in tip_fams or any(f == "other" for f in tip_fams):
        # discovery + unlabeled candidates share the green / grey lanes used for recovered tips
        if "discovery" in tip_fams:
            legend.append(Patch(color="#27AE60", label="Discovery candidate"))
        if "other" in tip_fams and "discovery" not in tip_fams:
            legend.append(Patch(color="#7F8C8D", label="Other / unlabeled"))
    if not legend:
        legend = [Patch(color="#7F8C8D", label="Tips")]
    ax.legend(handles=legend, loc="upper left", fontsize=8)
    _save(fig, out)


def plot_phylo_distance_heatmap(tree_path: Path, out: Path) -> None:
    if not tree_path.exists():
        return
    tree = Phylo.read(tree_path, "newick")
    tips = [t.name for t in tree.get_terminals()]
    n = len(tips)
    if n < 2 or n > 40:
        return
    dist = np.zeros((n, n))
    for i, a in enumerate(tips):
        for j, b in enumerate(tips):
            dist[i, j] = tree.distance(a, b)
    fig, ax = plt.subplots(figsize=(12, 10))
    im = ax.imshow(dist, cmap="viridis")
    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels([t[:20] for t in tips], rotation=90, fontsize=5)
    ax.set_yticklabels([t[:20] for t in tips], fontsize=5)
    ax.set_title("Phylogenetic distance matrix (gap junction panel)")
    plt.colorbar(im, ax=ax, label="Patristic distance")
    _save(fig, out)


# ---------------------------------------------------------------------------
# Synteny summary
# ---------------------------------------------------------------------------


def plot_synteny_job_summary(jobs: pd.DataFrame, out: Path) -> None:
    sub = jobs[jobs["status"].notna()].copy()
    fig, ax = plt.subplots(figsize=(9, 5))
    colors = ["#27AE60" if s == "complete" else "#E74C3C" for s in sub["status"]]
    ax.barh(sub["gene"], sub["synteny_block_plots"].fillna(0), color=colors, alpha=0.85)
    ax.set_xlabel("Synteny block plots generated")
    ax.set_title("SynVoy case study & reference panel jobs")
    _save(fig, out)


def plot_synteny_annotations(jobs: pd.DataFrame, out: Path) -> None:
    sub = jobs[jobs["status"] == "complete"].copy()
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    x = np.arange(len(sub))
    ax.bar(x - 0.2, sub["total_goi_annotations"].fillna(0), 0.4, label="Total GOI annotations", color="#3498DB")
    ax.bar(x + 0.2, sub["confident_goi_annotations"].fillna(0), 0.4, label="Confident", color="#27AE60")
    ax.set_xticks(x)
    ax.set_xticklabels(sub["gene"], rotation=25, ha="right")
    ax.set_title("SynVoy ortholog annotation yield by query gene")
    ax.legend()
    _save(fig, out)


# ---------------------------------------------------------------------------
# Extended comparisons (round 2)
# ---------------------------------------------------------------------------


def plot_case3_cross_clade(out: Path) -> None:
    csv_path = RESULTS_DIR / "synteny_case_studies" / "case3_cross_clade_inx2" / "case3_cross_clade_comparison.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    fig, axes = plt.subplots(1, 3, figsize=(12, 4))
    metrics = [
        ("confident_ortholog_species", "Confident orthologs"),
        ("probable_ortholog_species", "Probable orthologs"),
        ("total_goi_annotations", "GOI annotations"),
    ]
    for ax, (col, label) in zip(axes, metrics):
        ax.bar(df["clade"], df[col], color=[CLADE_COLORS.get(c, "#7F8C8D") for c in df["clade"]], alpha=0.85)
        ax.set_title(label)
        ax.set_ylabel("Count")
    fig.suptitle("Case 3: Inx2 cross-clade synteny (insect vs nematode)", y=1.02)
    _save(fig, out)


def plot_case2_ortholog_support(out: Path) -> None:
    csv_path = RESULTS_DIR / "synteny_case_studies" / "case2_ortholog_synteny" / "case2_ortholog_support.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    color_map = {
        "confident_ortholog": "#27AE60",
        "probable_ortholog": "#F39C12",
        "weak_goi_signal": "#E67E22",
        "flanking_synteny_only": "#95A5A6",
    }
    fig, ax = plt.subplots(figsize=(10, 5))
    labels = df["query_gene"] + " → " + df["target_species"].str.replace("Drosophila ", "D. ")
    colors = [color_map.get(s, "#BDC3C7") for s in df["ortholog_support"]]
    ax.barh(range(len(df)), df["best_region_score"], color=colors, alpha=0.85)
    ax.set_yticks(range(len(df)))
    ax.set_yticklabels(labels, fontsize=7)
    ax.set_xlabel("Best synteny region score")
    ax.set_title("Case 2: Inx2 vs shakB ortholog synteny support")
    legend = [Patch(color=c, label=k.replace("_", " ")) for k, c in color_map.items()]
    ax.legend(handles=legend, fontsize=7, loc="lower right")
    _save(fig, out)


def plot_inx2_ortholog_panel(gene_summary: pd.DataFrame, ref_inx: pd.DataFrame, out: Path) -> None:
    genes = [
        ("Drosophila melanogaster", "Inx2", "Q9V427"),
        ("Caenorhabditis elegans", "inx-2", "Q9U3K5"),
        ("Schistocerca americana", "inx2", "Q9XYN1"),
        ("Aedes aegypti", "shakB", "Q1DH70"),
        ("Anopheles gambiae", "shakB", "Q7PXN1"),
    ]
    rows = []
    for org, sym, acc in genes:
        g = gene_summary[(gene_summary["organism"] == org) & (gene_summary["reference_gene_symbol"].str.lower() == sym.lower())]
        r = ref_inx[(ref_inx["organism"] == org) & (ref_inx["uniprot_accession"] == acc)]
        if g.empty or r.empty:
            continue
        rows.append({
            "label": f"{org.split()[0][0]}. {sym}",
            "protein_length": r.iloc[0]["protein_length"],
            "exons": g.iloc[0]["mean_exon_count"],
            "gene_span_kb": g.iloc[0]["gene_span"] / 1000,
        })
    if not rows:
        return
    df = pd.DataFrame(rows)
    fig, axes = plt.subplots(1, 3, figsize=(13, 4))
    for ax, col, title in zip(axes, ["protein_length", "exons", "gene_span_kb"], ["Protein (aa)", "Exons", "Genomic span (kb)"]):
        ax.bar(df["label"], df[col], color="#2E86AB", alpha=0.85)
        ax.set_title(title)
        plt.setp(ax.xaxis.get_majorticklabels(), rotation=25, ha="right", fontsize=8)
    fig.suptitle("Inx2 / shakB ortholog panel across clades", y=1.02)
    _save(fig, out)


def plot_gjb2_bird_human(gene_summary: pd.DataFrame, ref_cnx: pd.DataFrame, out: Path) -> None:
    targets = ["Gallus gallus", "Homo sapiens", "Mus musculus"]
    sub = gene_summary[(gene_summary["reference_gene_symbol"] == "GJB2") & (gene_summary["organism"].isin(targets))]
    lens = ref_cnx[(ref_cnx["gene_symbol"] == "GJB2") & (ref_cnx["organism"].isin(targets))]
    merged = sub.merge(lens[["organism", "protein_length"]], on="organism")
    if merged.empty:
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    colors = ["#F1C40F", "#E74C3C", "#3498DB"]
    ax.bar(merged["organism"].str.replace("Gallus gallus", "Bird").str.replace("Homo sapiens", "Human").str.replace("Mus musculus", "Mouse"),
           merged["protein_length"], color=colors[: len(merged)], alpha=0.85)
    ax.set_ylabel("Protein length (aa)")
    ax.set_title("GJB2 (connexin 26): bird vs human vs mouse")
    _save(fig, out)


def plot_gja1_ortholog_trajectory(ref_cnx: pd.DataFrame, out: Path) -> None:
    sub = ref_cnx[ref_cnx["gene_symbol"] == "GJA1"].copy()
    sub["clade"] = sub["organism"].map(infer_clade)
    order = ["fish", "amphibian", "reptile", "bird", "mammal", "other"]
    sub["clade"] = pd.Categorical(sub["clade"], categories=order, ordered=True)
    sub = sub.sort_values("clade")
    fig, ax = plt.subplots(figsize=(10, 5))
    ax.plot(range(len(sub)), sub["protein_length"], "o-", color="#E67E22", linewidth=2, markersize=8)
    for i, row in sub.iterrows():
        idx = list(sub.index).index(i)
        ax.annotate(row["organism"].split()[0][0] + ". " + row["organism"].split()[1][:6], (idx, row["protein_length"]), fontsize=6, ha="center", va="bottom")
    ax.set_ylabel("GJA1 protein length (aa)")
    ax.set_xlabel("Species (clade-ordered)")
    ax.set_title("GJA1 ortholog length across vertebrate panel")
    _save(fig, out)


def plot_reference_vs_discovery(candidates: pd.DataFrame, ref_all: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5), sharey=True)
    for ax, family in zip(axes, ("innexin", "connexin")):
        ref = ref_all[ref_all["family"] == family]["protein_length"]
        disc = candidates[(candidates["family"] == family) & (candidates["rank_category"] != "rejected_false_positive")]["protein_length"]
        ax.hist(ref, bins=15, alpha=0.6, color=FAMILY_COLORS[family], label=f"Reference (n={len(ref)})", density=True)
        ax.hist(disc, bins=15, alpha=0.5, color="#27AE60", label=f"Discovery accepted (n={len(disc)})", density=True)
        ax.set_title(family.capitalize())
        ax.set_xlabel("Protein length (aa)")
        ax.legend(fontsize=7)
    axes[0].set_ylabel("Density")
    fig.suptitle("Reference vs discovery candidate lengths", y=1.02)
    _save(fig, out)


def plot_tm_helix_reference(structure: pd.DataFrame, ref_all: pd.DataFrame, out: Path) -> None:
    tm_counts = structure[structure["track"] == "transmembrane"].groupby("uniprot_accession").size()
    merged = ref_all.copy()
    merged["tm_helices"] = merged["uniprot_accession"].map(tm_counts).fillna(0).astype(int)
    merged = merged[merged["tm_helices"] > 0]
    if merged.empty:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    for family in ("innexin", "connexin"):
        sub = merged[merged["family"] == family]
        ax.scatter(sub["protein_length"], sub["tm_helices"], c=FAMILY_COLORS[family], alpha=0.7, s=50, label=family)
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("Predicted TM helices")
    ax.set_title("TM topology vs protein length (annotated subset)")
    ax.legend()
    _save(fig, out)


def plot_exon_length_distribution(out: Path) -> None:
    exon_path = EXON_DIR / "exon_coordinates.csv"
    if not exon_path.exists():
        return
    df = pd.read_csv(exon_path)
    fig, ax = plt.subplots(figsize=(8, 5))
    for family in ("innexin", "connexin"):
        sub = df[df["family"] == family]["exon_length"]
        ax.hist(sub, bins=30, alpha=0.5, color=FAMILY_COLORS[family], label=family, density=True)
    ax.set_xlabel("Exon length (bp)")
    ax.set_ylabel("Density")
    ax.set_title("Coding exon length distribution by family")
    ax.legend()
    _save(fig, out)


def plot_transcript_complexity(gene_summary: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for family in ("innexin", "connexin"):
        sub = gene_summary[gene_summary["family"] == family]
        ax.scatter(sub["transcript_count"], sub["mean_exon_count"], c=FAMILY_COLORS[family], alpha=0.6, s=45, label=family)
    ax.set_xlabel("Transcript isoforms per gene")
    ax.set_ylabel("Mean exon count")
    ax.set_title("Transcript complexity vs exon architecture")
    ax.legend()
    _save(fig, out)


def plot_rank_category_breakdown(candidates: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(9, 5))
    cats = ["high_confidence", "possible", "weak", "rejected"]
    x = np.arange(2)
    width = 0.18
    for i, cat in enumerate(cats):
        vals = []
        for family in ("innexin", "connexin"):
            sub = candidates[candidates["family"] == family]
            if cat == "rejected":
                n = len(sub[sub["rank_category"] == "rejected_false_positive"])
            elif cat == "high_confidence":
                n = len(sub[sub["rank_category"].str.contains("high_confidence", na=False)])
            elif cat == "weak":
                n = len(sub[sub["rank_category"] == "weak_manual_review"])
            else:
                n = len(sub[sub["rank_category"].str.contains("possible", na=False)])
            vals.append(n)
        ax.bar(x + i * width, vals, width, label=cat.replace("_", " "))
    ax.set_xticks(x + 1.5 * width)
    ax.set_xticklabels(["Innexin", "Connexin"])
    ax.set_ylabel("Candidates")
    ax.set_title("Discovery rank category breakdown")
    ax.legend(fontsize=7)
    _save(fig, out)


def plot_nematode_innexin_panel(ref_inx: pd.DataFrame, out: Path) -> None:
    sub = ref_inx[ref_inx["organism"].str.contains("Caenorhabditis")].sort_values(["organism", "protein_length"])
    fig, ax = plt.subplots(figsize=(9, 6))
    labels = sub["organism"].str.replace("Caenorhabditis ", "C. ") + " | " + sub["gene_symbol"]
    ax.barh(labels, sub["protein_length"], color=["#9B59B6", "#8E44AD"] * (len(sub) // 2 + 1), alpha=0.85)
    ax.set_xlabel("Protein length (aa)")
    ax.set_title("Nematode innexin family — C. elegans vs C. briggsae")
    _save(fig, out)


def plot_synteny_goi_heatmap(out: Path) -> None:
    csv_path = METADATA_DIR / "synteny_genome_comparison.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    pivot = df.pivot_table(index="gene", columns="target_genome", values="goi_annotations", aggfunc="max").fillna(0)
    if pivot.empty:
        return
    fig, ax = plt.subplots(figsize=(12, 4))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlGnBu")
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels([c[:15] for c in pivot.columns], rotation=45, ha="right", fontsize=6)
    ax.set_title("GOI annotations per target genome (SynVoy heatmap)")
    plt.colorbar(im, ax=ax)
    _save(fig, out)


def plot_assembly_quality(out: Path) -> None:
    qpath = METADATA_DIR / "gap_junction_quality_table.csv"
    if not qpath.exists():
        return
    df = pd.read_csv(qpath)
    fig, axes = plt.subplots(1, 2, figsize=(11, 4))
    for ax, family in zip(axes, ("innexin", "connexin")):
        sub = df[df["family"] == family]
        ann = sub["annotation_available"].value_counts()
        ax.pie(ann.values, labels=[f"GFF: {k}" for k in ann.index], autopct="%1.0f%%", colors=["#27AE60", "#E74C3C"][: len(ann)])
        ax.set_title(f"{family.capitalize()} reference assembly annotation")
    fig.suptitle("Genome annotation availability (reference panel)", y=1.02)
    _save(fig, out)


def plot_connexin_batch_discovery(validation: pd.DataFrame, out: Path) -> None:
    batches = validation[validation["family"] == "connexin"]["discovery_batch"].unique()
    if len(batches) < 2:
        sub = validation[validation["family"] == "connexin"]
        fig, ax = plt.subplots(figsize=(8, 5))
        ax.barh(sub["organism"].str[:20], sub["accepted_total"], color="#E67E22", alpha=0.8)
        ax.set_xlabel("Accepted candidates")
        ax.set_title("Connexin discovery acceptance by species")
        _save(fig, out)
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    agg = validation[validation["family"] == "connexin"].groupby("discovery_batch")[["candidate_total", "accepted_total", "high_confidence"]].sum()
    x = np.arange(len(agg))
    ax.bar(x - 0.25, agg["candidate_total"], 0.25, label="Total", color="#BDC3C7")
    ax.bar(x, agg["accepted_total"], 0.25, label="Accepted", color="#E67E22")
    ax.bar(x + 0.25, agg["high_confidence"], 0.25, label="High confidence", color="#27AE60")
    ax.set_xticks(x)
    ax.set_xticklabels([b.replace("connexin_", "").replace("_discovery", "") for b in agg.index], rotation=15, ha="right")
    ax.set_title("Connexin discovery: batch comparison")
    ax.legend()
    _save(fig, out)


def plot_innexin_clade_discovery(thesis: pd.DataFrame, out: Path) -> None:
    sub = thesis[(thesis["family"] == "innexin") & (thesis["panel_role"].str.contains("discovery"))]
    if sub.empty:
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    sub = sub.sort_values("discovery_accepted", ascending=True)
    colors = ["#27AE60" if v > 0 else "#E74C3C" for v in sub["discovery_accepted"].fillna(0)]
    ax.barh(sub["organism"].str[:25], sub["discovery_accepted"].fillna(0), color=colors, alpha=0.85)
    ax.set_xlabel("Accepted innexin candidates")
    ax.set_title("Innexin discovery success by target species")
    _save(fig, out)


def plot_shakB_insect_comparison(ref_inx: pd.DataFrame, out: Path) -> None:
    sub = ref_inx[ref_inx["gene_symbol"].str.lower() == "shakb"].sort_values("protein_length")
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.bar(sub["organism"].str.replace("Drosophila melanogaster", "D. melanogaster").str.replace("Aedes aegypti", "A. aegypti").str.replace("Anopheles gambiae", "A. gambiae"),
           sub["protein_length"], color="#1ABC9C", alpha=0.85)
    ax.set_ylabel("Protein length (aa)")
    ax.set_title("shakB innexin orthologs — insect panel")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=15, ha="right")
    _save(fig, out)


def plot_flanking_gene_counts(out: Path) -> None:
    csv_path = METADATA_DIR / "synteny_home_flanking_genes.csv"
    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    fig, ax = plt.subplots(figsize=(8, 5))
    for gene, grp in df.groupby("gene"):
        ax.scatter(grp["flanking_gene_count"], grp["goi_count"], s=60, alpha=0.7, label=gene)
    ax.set_xlabel("Flanking genes in synteny block")
    ax.set_ylabel("GOI copies in block")
    ax.set_title("Home synteny block: flanking genes vs GOI count")
    ax.legend(fontsize=8)
    _save(fig, out)


def plot_secondary_structure_fraction(structure: pd.DataFrame, ref_all: pd.DataFrame, out: Path) -> None:
    rows = []
    for acc, grp in structure[structure["track"] == "secondary_structure"].groupby("uniprot_accession"):
        plen = ref_all[ref_all["uniprot_accession"] == acc]["protein_length"]
        if plen.empty:
            continue
        length = int(plen.iloc[0])
        helix = sum(int(r["aa_end"]) - int(r["aa_start"]) + 1 for _, r in grp[grp["feature_label"] == "Helix"].iterrows())
        sheet = sum(int(r["aa_end"]) - int(r["aa_start"]) + 1 for _, r in grp[grp["feature_label"] == "Sheet"].iterrows())
        fam = ref_all[ref_all["uniprot_accession"] == acc]["family"].iloc[0]
        rows.append({"accession": acc, "family": fam, "helix_frac": helix / length, "sheet_frac": sheet / length})
    if not rows:
        return
    df = pd.DataFrame(rows)
    fig, ax = plt.subplots(figsize=(8, 5))
    for family in ("innexin", "connexin"):
        sub = df[df["family"] == family]
        ax.scatter(sub["helix_frac"], sub["sheet_frac"], c=FAMILY_COLORS[family], alpha=0.7, s=50, label=family)
    ax.set_xlabel("Helix fraction")
    ax.set_ylabel("Sheet fraction")
    ax.set_title("Secondary structure composition (ProtT5)")
    ax.legend()
    _save(fig, out)


def plot_inx2_vs_inx2_nematode_insect(gene_summary: pd.DataFrame, ref_inx: pd.DataFrame, out: Path) -> None:
    pairs = [
        ("Drosophila melanogaster", "Inx2"),
        ("Caenorhabditis elegans", "inx-2"),
    ]
    metrics = []
    for org, sym in pairs:
        g = gene_summary[(gene_summary["organism"] == org) & (gene_summary["reference_gene_symbol"].str.lower() == sym.lower())]
        r = ref_inx[(ref_inx["organism"] == org) & (ref_inx["gene_symbol"].str.lower() == sym.lower())]
        if g.empty or r.empty:
            continue
        metrics.append({
            "group": "Insect (Dmel)" if "Drosophila" in org else "Nematode (Cele)",
            "protein_length": r.iloc[0]["protein_length"],
            "exons": g.iloc[0]["mean_exon_count"],
            "span_kb": g.iloc[0]["gene_span"] / 1000,
        })
    if len(metrics) < 2:
        return
    df = pd.DataFrame(metrics)
    fig, axes = plt.subplots(1, 3, figsize=(10, 3.5))
    colors = ["#1ABC9C", "#9B59B6"]
    for ax, col, ylab in zip(axes, ["protein_length", "exons", "span_kb"], ["Protein (aa)", "Exons", "Span (kb)"]):
        ax.bar(df["group"], df[col], color=colors, alpha=0.85)
        ax.set_title(ylab)
    fig.suptitle("Inx2: insect vs nematode architecture (Case 3 anchor)", y=1.05)
    _save(fig, out)


# ---------------------------------------------------------------------------
# HTML gallery (presentation-ready)
# ---------------------------------------------------------------------------

GALLERY_SECTIONS: tuple[tuple[str, str, str], ...] = (
    ("overview", "Overview", "Scope of the species panel and pipeline status."),
    ("discovery", "Discovery", "Genome-wide candidate search and validation outcomes."),
    ("comparison", "Innexin vs connexin", "Cross-family architecture and discovery metrics."),
    ("phylogeny", "Phylogeny", "Evolutionary relationships among gap-junction proteins."),
    ("exon_maps", "Exon architecture", "Genomic and compressed exon–intron maps (±180 kb neighbour windows)."),
    ("domains", "Domains & structure", "Pfam, TM topology, disorder, and secondary structure."),
    ("protein_features", "Protein feature maps", "Domain tracks for key orthologs."),
    ("innexin", "Innexin focus", "Invertebrate gap-junction gene family in detail."),
    ("connexin", "Connexin focus", "Vertebrate connexin landscape and ortholog trends."),
    ("synteny", "Synteny (SynVoy)", "Micro-synteny around innexin loci across genomes."),
    ("cross_clade", "Cross-clade cases", "Insect vs nematode architecture comparisons."),
    ("ortholog", "Ortholog support", "Case-study ortholog evidence (Inx2, shakB)."),
)

# Selected findings for the presenter page.
# Each entry: rel path, headline, caption, optional speaker note.
KEY_FINDINGS: tuple[dict[str, str], ...] = (
    {
        "rel": "imported/exon_maps/innexin_genomic_exon_map.png",
        "headline": "Innexin exon counts differ a lot across species",
        "caption": (
            "C. elegans innexins reach up to 9 exons; Drosophila shakB has only 1 exon; "
            "ogre spans 6–8 exons. Same protein family, very different gene structures."
        ),
        "script": (
            "Innexin gene "
            "architecture is not conserved. Nematode inx genes are multi-exon, Drosophila shakB "
            "is a single exon, and ogre is huge. Architecture alone cannot define orthology."
        ),
    },
    {
        "rel": "innexin/04_dmel_paralog_lengths.png",
        "headline": "Within Drosophila, ogre is 3× larger than other innexin paralogs",
        "caption": (
            "ogre gene span ~17 kb vs Inx2/Inx7 ~2–5 kb. Size heterogeneity exists even "
            "among paralogs in one species."
        ),
        "script": (
            "Even within D. melanogaster, innexin paralogs differ dramatically. ogre is "
            "the giant — long gene, many exons — while Inx2 and Inx7 are compact. "
            "That difference motivates the synteny cluster comparison."
        ),
    },
    {
        "rel": "synteny/04_case1_cluster_sharing.png",
        "headline": "Dmel innexin paralogs share 100% of flanking genes",
        "caption": (
            "Case 1: ogre, Inx2, and Inx3 at both loci share all 25–27 neighbouring genes "
            "with zero unique flanking genes — classic tandem duplication cluster."
        ),
        "script": (
            "At both the proximal and distal "
            "clusters, ogre, Inx2, and Inx3 share every flanking gene — 25 to 27 genes, "
            "zero unique neighbours. Consistent with local tandem duplication, "
            "not independent insertions."
        ),
    },
    {
        "rel": "ortholog/01_case2_support.png",
        "headline": "Inx2 has probable orthologs in 3 Drosophila species; shakB only 1 confident hit",
        "caption": (
            "Case 2 SynVoy: Inx2 → probable ortholog in D. limbata, D. helvetica, D. busckii. "
            "shakB → confident ortholog only in D. busckii; flanking-only in D. helvetica."
        ),
        "script": (
            "Cross-species orthology is gene-specific. Inx2 gets probable support in "
            "three Drosophila species. shakB is much stricter — only D. busckii has a confident "
            "ortholog; in D. helvetica synteny is present but no clear GOI match."
        ),
    },
    {
        "rel": "ortholog/03_shakB_insects.png",
        "headline": "shakB length and architecture differ across insects",
        "caption": (
            "shakB orthologs in Aedes, Anopheles, and Drosophila differ in protein length "
            "and exon count — not a simple copy-paste ortholog."
        ),
        "script": (
            "Zooming into shakB across insects: protein length and exon structure are "
            "not uniform. Anopheles shakB has up to 8 exons, Aedes has 10 — very different "
            "from Drosophila's single-exon shakB."
        ),
    },
    {
        "rel": "cross_clade/01_case3_comparison.png",
        "headline": "Insect Inx2: 58 synteny hits vs nematode inx-2: only 6",
        "caption": (
            "Case 3: SynVoy across 5 insect genomes yields 58 GOI annotations for Inx2; "
            "nematode inx-2 panel yields only 6 — synteny signal is phylum-specific."
        ),
        "script": (
            "Cross-clade comparison: insect Inx2 produces 58 gene-of-interest annotations "
            "across five genomes, but C. elegans inx-2 only six. Synteny conservation does not "
            "transfer across phyla — we cannot use insect synteny rules for nematodes."
        ),
    },
    {
        "rel": "cross_clade/02_inx2_architecture.png",
        "headline": "Inx2 vs inx-2: same name, different exon architecture",
        "caption": (
            "Drosophila Inx2: 1–3 exons, ~5480 bp gene span. C. elegans inx-2: 4 exons, "
            "1727 bp. Ortholog naming does not imply structural similarity."
        ),
        "script": (
            "The architecture bar chart makes this concrete: Drosophila Inx2 is 2-exon "
            "and compact, C. elegans inx-2 has 4 exons. Calling both Inx2 is historical naming — "
            "the gene structures diverged."
        ),
    },
    {
        "rel": "imported/protein_features/shakB_protein_feature_map.png",
        "headline": "shakB keeps canonical innexin TM topology across a multi-exon, multi-transcript locus",
        "caption": (
            "Predicted TM helices and disorder along shakB. The locus is intron-rich and splice-capable "
            "(many annotated transcripts in mosquito/fly panels), yet still encodes a plausible gap-junction channel protein."
        ),
        "script": (
            "shakB is a classic multi-isoform innexin locus — coding-region introns, not a single-exon gene. "
            "The protein feature map still shows the expected transmembrane topology. Exon count and domain architecture "
            "are related but not the same thing."
        ),
    },
    {
        "rel": "discovery/02_innexin_candidates.png",
        "headline": (
            "True rescues: 3 species where discovery returned zero accepted "
            "(≥28 loci, floor); plus first-pass curator hits elsewhere"
        ),
        "caption": (
            "Rescue means prior discovery ran and accepted=0 — Rotaria macrura, "
            "Brachionus manjavacas, Abra alba — not every species with prior_discovery=not_run. "
            "Most curator positives are first-time probes of genomes never run through discovery. "
            "Locus counts are floors (Dmel: miniprot-only 3/8 → region locate 7–8/8). "
            "Counts near the query-pack size can track queries; region locate can exceed the pack "
            "(e.g. Acanthocardia ≥22 with 15 queries)."
        ),
        "script": (
            "The first automated innexin search looked sparse. Correct framing: only three species "
            "are true rescues (discovery ran, accepted=0). The rest of the curator yield is "
            "first-pass probing of not_run genomes, plus floors — not a full paralog census."
        ),
    },
    {
        "rel": "comparison/08_intron_architecture.png",
        "headline": "Connexins have longer introns than innexins at similar exon counts",
        "caption": (
            "Scatter of exon count vs mean intron length (log scale). Connexin genes occupy "
            "a different intron-size regime than innexins."
        ),
        "script": (
            "At comparable exon numbers, connexin introns are systematically longer. "
            "This supports treating the two families as separate architectural lineages, not "
            "just different species panels."
        ),
    },
    {
        "rel": "connexin/06_gja1_trajectory.png",
        "headline": "GJA1 protein length conserved across vertebrates (~380 aa)",
        "caption": (
            "GJA1 ortholog trajectory across mammals, birds, fish. Connexin 43 is highly "
            "conserved compared to other connexin paralogs."
        ),
        "script": (
            "On the connexin side, GJA1 is remarkably conserved — about 380 amino acids "
            "across vertebrates. That contrasts with the architectural diversity among innexins."
        ),
    },
)

FINDING_RELS: frozenset[str] = frozenset(f["rel"] for f in KEY_FINDINGS)

# Full presenter gallery (all sections, captions, scripts — unchanged layout).
PRESENTER_KEY_FIGURES: frozenset[str] = frozenset(
    {
        "overview/01_panel_composition.png",
        "overview/02_discovery_funnel.png",
        "comparison/01_protein_length_boxplot.png",
        "phylogeny/01_tree_colored.png",
        "phylogeny/03_innexin_reference_tree.png",
        "synteny/03_refpanel_status.png",
        "synteny/04_case1_cluster_sharing.png",
        "imported/exon_maps/innexin_genomic_exon_map.png",
        "imported/protein_features/shakB_protein_feature_map.png",
        "cross_clade/02_inx2_architecture.png",
        "discovery/01_validation_by_species.png",
    }
)

PRESENTER_TALK_TRACK: tuple[tuple[str, str, str], ...] = (
    ("Overview", "#overview", "Introduce species panel and annotation gaps"),
    ("Discovery", "#discovery", "Show candidate funnel and per-species validation"),
    ("Comparison", "#comparison", "Innexin vs connexin — core thesis contrast"),
    ("Phylogeny", "#phylogeny", "Separate evolutionary histories of both families"),
    ("Exon maps", "#exon_maps", "Gene architecture slides (±180 kb neighbour windows)"),
    ("Synteny", "#synteny", "SynVoy reference panels + Dmel cluster case"),
)

PRESENTER_FIGURE_CAPTIONS: dict[str, str] = {
    "overview/01_panel_composition.png": (
        "Which species are in the thesis panel and whether they serve as curated references "
        "or discovery-only genomes. Sets the taxonomic scope (insects, nematodes, vertebrates)."
    ),
    "overview/02_discovery_funnel.png": (
        "How many raw candidates survive each filtering step (length, transmembrane helices, "
        "reference similarity). Shows that discovery is intentionally strict."
    ),
    "overview/04_annotation_status.png": (
        "Which genomes already have gap-junction annotation in NCBI vs. species where we "
        "must discover genes de novo — motivates the custom pipeline."
    ),
    "comparison/01_protein_length_boxplot.png": (
        "Core biological comparison: innexins (invertebrate) vs connexins (vertebrate) protein "
        "lengths. Connexins tend to be longer and more variable in mammals."
    ),
    "comparison/06_acceptance_rates.png": (
        "Accepted vs rejected discovery candidates by family. Most connexin hits come from "
        "large evolutionary discovery batches; innexin discovery is harder in poorly annotated genomes."
    ),
    "discovery/01_validation_by_species.png": (
        "Per-species yield of validated candidates. Highlights which discovery genomes produced "
        "usable gap-junction hits (and which did not, e.g. Rotaria)."
    ),
    "discovery/02_innexin_candidates.png": (
        "All non-rejected innexin candidates with protein length and confidence colour. "
        "Use to discuss manual review priorities."
    ),
    "phylogeny/01_tree_colored.png": (
        "Combined phylogeny of reference innexins and connexins. The two families form "
        "separate clades — they are homologous channel proteins but not interchangeable orthologs."
    ),
    "phylogeny/03_innexin_reference_tree.png": (
        "Innexin-only reference panel used for SynVoy synteny (Drosophila, C. elegans, etc.). "
        "Shows which paralogs anchor cross-species comparisons."
    ),
    "synteny/03_refpanel_status.png": (
        "Status of the three completed SynVoy reference-panel jobs (Inx2, shakB, inx-2). "
        "Key milestone: all reference panels finished with synteny plots."
    ),
    "synteny/04_case1_cluster_sharing.png": (
        "Case study: do Drosophila innexin paralogs share flanking genes? Shared synteny "
        "suggests local duplication; unique flanking genes suggest relocation or divergence."
    ),
    "cross_clade/02_inx2_architecture.png": (
        "Case 3: Inx2 exon count and protein length in insects vs nematodes. Tests whether "
        "the same ortholog name implies comparable gene architecture across phyla."
    ),
    "ortholog/02_inx2_panel.png": (
        "Inx2 ortholog panel across the reference species — length and exon structure side by side."
    ),
    "imported/exon_maps/innexin_genomic_exon_map.png": (
        "Genomic-scale exon map for innexins: each row is a species/gene, boxes are exons. "
        "Directly comparable to classic FOXP2 presentation slides."
    ),
    "imported/exon_maps/connexin_genomic_exon_map_part1.png": (
        "Connexin exon map (part 1): vertebrate connexins often have more exons and longer "
        "introns than innexins."
    ),
    "imported/protein_features/shakB_protein_feature_map.png": (
        "Domain architecture of Drosophila shakB: transmembrane segments, disorder, and motifs "
        "along the protein sequence."
    ),
    "imported/protein_features/GJA1_human_protein_feature_map.png": (
        "Human GJA1 (connexin 43) feature map for comparison with invertebrate innexins."
    ),
    "domains/domain_mv_pfam_connexin.png": (
        "Monument-valley plot: Pfam Connexin domain coverage across the reference panel."
    ),
    "domains/domain_mv_tm_topology.png": (
        "Predicted transmembrane topology coverage — gap-junction proteins should show "
        "4–6 TM helices."
    ),
    "innexin/04_dmel_paralog_lengths.png": (
        "All D. melanogaster innexin paralog lengths — basis for the Dmel synteny cluster case study."
    ),
}
for _finding in KEY_FINDINGS:
    PRESENTER_FIGURE_CAPTIONS[_finding["rel"]] = _finding["caption"]

PRESENTER_SPEAKER_SCRIPT: dict[str, str] = {
    "overview/01_panel_composition.png": (
        "The panel covers gap-junction proteins across 134 species — insects, nematodes, "
        "and vertebrates. Some species have curated reference genes; others need de novo discovery "
        "because NCBI lacks annotation."
    ),
    "overview/02_discovery_funnel.png": (
        "We search genomes with MMseqs/miniprot, then filter by protein length, transmembrane "
        "helices, and similarity to known innexins/connexins. Most hits are rejected — the pipeline "
        "is intentionally strict."
    ),
    "comparison/01_protein_length_boxplot.png": (
        "Innexins and connexins are homologous channel families but not orthologs. Connexins "
        "are generally longer, especially in mammals — so we analyse them separately."
    ),
    "phylogeny/01_tree_colored.png": (
        "The phylogeny confirms two distinct clades. Any cross-family comparison is structural, "
        "not evolutionary orthology."
    ),
    "phylogeny/03_innexin_reference_tree.png": (
        "These innexins anchor our SynVoy synteny panels — Inx2 and shakB in Drosophila, "
        "inx-2 in C. elegans."
    ),
    "synteny/03_refpanel_status.png": (
        "All three reference-panel SynVoy jobs finished successfully — this was a major "
        "computational milestone (OOM fixes, cached genomes, laptop_safe profile)."
    ),
    "synteny/04_case1_cluster_sharing.png": (
        "For Drosophila innexin paralogs we ask: do they sit in duplicated blocks with shared "
        "neighbours? Shared flanking genes support local duplication."
    ),
    "imported/exon_maps/innexin_genomic_exon_map.png": (
        "Exon maps follow the FOXP2 presentation style — each row is a gene, boxes are exons. "
        "Innexins vary in exon count but share the core TM-domain architecture."
    ),
    "imported/protein_features/shakB_protein_feature_map.png": (
        "Protein feature tracks show TM helices and disordered regions along shakB — our main "
        "insect electrical-synapse innexin example."
    ),
    "cross_clade/02_inx2_architecture.png": (
        "Inx2 in insects vs nematodes — same gene name, different phyla. We check whether "
        "exon count and protein length are comparable before calling them orthologs."
    ),
    "discovery/01_validation_by_species.png": (
        "Discovery yield differs by species. Curator re-annotation rescued several "
        "zero-accepted cases — notably Rotaria macrura and Brachionus manjavacas — and added "
        "new arthropod and mollusc innexin loci."
    ),
}
for _finding in KEY_FINDINGS:
    if _finding.get("script"):
        PRESENTER_SPEAKER_SCRIPT[_finding["rel"]] = _finding["script"]


def _gallery_stats() -> dict[str, str]:
    stats = {
        "species": "134",
        "candidates": "359",
        "accepted": "144",
        "synvoy": "3/3",
        "curator_present_spp": "21",
        "curator_present_loci": "250",
        # True rescue = prior discovery ran and accepted=0 (not the same as first-pass not_run).
        "curator_rescued_spp": "3",
        "curator_rescued_loci": "28",
        "curator_firstpass_spp": "14",
        "curator_firstpass_loci": "178",
        # Legacy combined bucket (rescued + first-pass); prefer rescued/firstpass in copy.
        "curator_new_spp": "17",
        "curator_new_loci": "206",
    }
    try:
        cand = pd.read_csv(METADATA_DIR / "gap_junction_candidates_master.csv")
        panel = pd.read_csv(METADATA_DIR / "species_panel_master.csv")
        stats["species"] = str(len(panel))
        stats["candidates"] = str(len(cand))
        stats["accepted"] = str(len(cand[cand["rank_category"] != "rejected_false_positive"]))
        summary = RESULTS_DIR / "gene_curator_probe" / "innexin_search_summary.csv"
        if summary.exists():
            cur = pd.read_csv(summary)
            present_spp = cur[cur["present_loci"].fillna(0).astype(int) > 0]
            stats["curator_present_spp"] = str(len(present_spp))
            stats["curator_present_loci"] = str(int(present_spp["present_loci"].sum()))

            def _prior_accepted(row: pd.Series) -> int | None:
                prior = str(row.get("prior_discovery") or "")
                if prior == "not_run":
                    return None
                try:
                    return int(row.get("prior_accepted"))
                except (TypeError, ValueError):
                    return None

            rescued = present_spp[present_spp.apply(lambda r: _prior_accepted(r) == 0, axis=1)]
            firstpass = present_spp[present_spp.apply(lambda r: _prior_accepted(r) is None, axis=1)]
            stats["curator_rescued_spp"] = str(len(rescued))
            stats["curator_rescued_loci"] = str(int(rescued["present_loci"].sum()) if len(rescued) else 0)
            stats["curator_firstpass_spp"] = str(len(firstpass))
            stats["curator_firstpass_loci"] = str(int(firstpass["present_loci"].sum()) if len(firstpass) else 0)
            # Combined "hit for the first time in this project" (rescued + first-pass)
            stats["curator_new_spp"] = str(len(rescued) + len(firstpass))
            stats["curator_new_loci"] = str(
                int(rescued["present_loci"].sum() if len(rescued) else 0)
                + int(firstpass["present_loci"].sum() if len(firstpass) else 0)
            )
    except (OSError, pd.errors.EmptyDataError, KeyError, ValueError):
        pass
    return stats


def _clean_plot_title(stem: str) -> str:
    """Remove vendor tokens from user-visible plot titles."""
    title = stem.replace("_", " ")
    title = re.sub(r"\brostlab\b", "", title, flags=re.IGNORECASE)
    return re.sub(r"\s+", " ", title).strip()


def _collect_imported_plots(out_dir: Path) -> list[tuple[str, str, Path]]:
    imported: list[tuple[str, str, Path]] = []
    exon_dir = out_dir / "imported" / "exon_maps"
    feat_dir = out_dir / "imported" / "protein_features"
    for path in sorted(exon_dir.glob("*.png")) if exon_dir.is_dir() else []:
        imported.append(("exon_maps", _clean_plot_title(path.stem), path))
    for path in sorted(feat_dir.glob("*.png")) if feat_dir.is_dir() else []:
        imported.append(("protein_features", _clean_plot_title(path.stem), path))
    return imported


def _escape_html(text: str) -> str:
    return text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;").replace('"', "&quot;")


def _prepare_gallery_sections(
    plot_records: list[tuple[str, str, Path]], out_dir: Path
) -> dict[str, list[tuple[str, str, Path]]]:
    sections: dict[str, list[tuple[str, str, Path]]] = defaultdict(list)
    for section, title, path in plot_records:
        if path.exists():
            sections[section].append((title, path))
    for section, title, path in _collect_imported_plots(out_dir):
        sections[section].append((title, path))
    return sections


_GALLERY_CSS = """
:root{--blue:#2E86AB;--orange:#E67E22;--green:#27AE60;--ink:#2C3E50;--muted:#5d6d7e;--bg:#f4f7fb;--card:#fff}
*{box-sizing:border-box}
body{font-family:'Segoe UI',system-ui,sans-serif;margin:0;background:var(--bg);color:var(--ink);line-height:1.55}
.layout{display:grid;grid-template-columns:240px 1fr;min-height:100vh}
nav{position:sticky;top:0;height:100vh;overflow:auto;padding:1.25rem 1rem;background:#1a2836;color:#ecf0f1;border-right:1px solid #2c3e50}
nav h2{font-size:.75rem;text-transform:uppercase;letter-spacing:.08em;color:#95a5a6;margin:0 0 .75rem}
nav a{display:block;color:#ecf0f1;text-decoration:none;padding:.35rem .5rem;border-radius:6px;font-size:.88rem;margin-bottom:.15rem}
nav a:hover{background:rgba(255,255,255,.08)}
nav .switch{margin-top:1.5rem;padding-top:1rem;border-top:1px solid #34495e;font-size:.8rem;color:#95a5a6}
nav .switch a{color:#5dade2;font-weight:600}
main{padding:2rem 2.5rem 4rem;max-width:1100px}
.hero{background:linear-gradient(135deg,var(--blue),var(--orange));color:#fff;padding:2rem 2.25rem;border-radius:16px;margin-bottom:2rem;box-shadow:0 8px 32px rgba(46,134,171,.25)}
.hero.presenter{background:linear-gradient(135deg,#1a5276,#117a65)}
.hero h1{margin:0 0 .5rem;font-size:1.85rem;font-weight:700}
.hero .subtitle{opacity:.95;font-size:1.05rem;margin:0 0 1.25rem}
.stats{display:flex;flex-wrap:wrap;gap:.75rem}
.stat{background:rgba(255,255,255,.18);backdrop-filter:blur(4px);padding:.5rem 1rem;border-radius:999px;font-size:.85rem}
.stat b{font-size:1.1rem;margin-right:.25rem}
.guide{background:#fff;border-left:4px solid var(--green);padding:1.25rem 1.5rem;border-radius:0 12px 12px 0;margin-bottom:2.5rem;box-shadow:0 2px 12px rgba(0,0,0,.06)}
.guide h2{margin:0 0 .75rem;font-size:1.15rem;color:var(--green)}
.guide ol{margin:0;padding-left:1.25rem}
.guide li{margin-bottom:.4rem;color:var(--muted)}
.guide a{color:var(--blue)}
section{margin-bottom:3rem;scroll-margin-top:1rem}
section > h2{font-size:1.45rem;color:var(--blue);border-bottom:2px solid var(--blue);padding-bottom:.4rem;margin:0 0 .5rem}
.section-intro{color:var(--muted);margin:0 0 1.25rem;font-size:.98rem;max-width:72ch}
.grid{display:grid;grid-template-columns:repeat(auto-fill,minmax(340px,1fr));gap:1.25rem}
.card{background:var(--card);border-radius:12px;padding:1rem 1.1rem 1.15rem;box-shadow:0 2px 14px rgba(0,0,0,.07);border:1px solid #e8edf3;transition:transform .15s,box-shadow .15s}
.card:hover{transform:translateY(-2px);box-shadow:0 6px 20px rgba(0,0,0,.1)}
.card.key{border-color:var(--green);box-shadow:0 2px 14px rgba(39,174,96,.15)}
.card img{width:100%;border-radius:8px;cursor:zoom-in;background:#f8f9fa}
.card h3{margin:.75rem 0 .35rem;font-size:1rem;color:var(--ink)}
.card h3.peer-title{text-align:center;font-weight:600;margin:.85rem 0 0;font-size:.92rem}
.card .caption{margin:.5rem 0 0;font-size:.88rem;color:var(--muted)}
.card .script{margin:.65rem 0 0;padding:.65rem .75rem;background:#eafaf1;border-radius:8px;font-size:.86rem;color:#1e5631;border-left:3px solid var(--green)}
.badge{display:inline-block;background:var(--green);color:#fff;font-size:.65rem;font-weight:700;text-transform:uppercase;letter-spacing:.05em;padding:.2rem .5rem;border-radius:4px;margin-bottom:.35rem}
.finding-num{display:inline-flex;align-items:center;justify-content:center;width:2rem;height:2rem;border-radius:50%;background:var(--blue);color:#fff;font-weight:700;font-size:.9rem;margin-bottom:.5rem}
.finding-headline{font-size:1.05rem;font-weight:600;color:var(--blue);margin:.5rem 0 .75rem;line-height:1.35}
.findings-list{display:flex;flex-direction:column;gap:2rem;max-width:920px}
.finding-card{background:var(--card);border-radius:14px;padding:1.25rem 1.5rem;box-shadow:0 4px 20px rgba(0,0,0,.08);border:1px solid #e8edf3}
.finding-card img{width:100%;border-radius:10px;cursor:zoom-in;background:#f8f9fa;margin-top:.5rem}
details.appendix{margin-top:3rem;background:#fff;border-radius:12px;padding:1rem 1.5rem;box-shadow:0 2px 12px rgba(0,0,0,.06)}
details.appendix summary{cursor:pointer;font-weight:600;color:var(--muted);padding:.5rem 0}
details.appendix summary:hover{color:var(--blue)}
.modal{display:none;position:fixed;inset:0;background:rgba(0,0,0,.88);z-index:1000;align-items:center;justify-content:center;padding:2rem;cursor:zoom-out}
.modal.open{display:flex}
.modal img{max-width:95vw;max-height:90vh;border-radius:8px;box-shadow:0 8px 40px rgba(0,0,0,.5)}
@media(max-width:900px){.layout{grid-template-columns:1fr}nav{position:relative;height:auto}main{padding:1.25rem}}
"""


def _path_lookup(
    sections: dict[str, list[tuple[str, str, Path]]], out_dir: Path
) -> dict[str, tuple[str, Path]]:
    lookup: dict[str, tuple[str, Path]] = {}
    for items in sections.values():
        for title, path in items:
            rel = path.relative_to(out_dir).as_posix()
            lookup[rel] = (title, path)
    return lookup


def _ordered_key_findings(
    sections: dict[str, list[tuple[str, str, Path]]],
    out_dir: Path,
) -> list[dict[str, str]]:
    """Return KEY_FINDINGS in the same order figures appear in index.html."""
    by_rel = {f["rel"]: f for f in KEY_FINDINGS}
    ordered: list[dict[str, str]] = []
    for sid, _label, _intro in GALLERY_SECTIONS:
        if sid not in sections:
            continue
        for _title, path in sections[sid]:
            rel = path.relative_to(out_dir).as_posix()
            if rel in by_rel:
                ordered.append(by_rel[rel])
    return ordered


def _render_peer_gallery_html(
    sections: dict[str, list[tuple[str, str, Path]]],
    out_dir: Path,
) -> str:
    stats = _gallery_stats()
    parts: list[str] = [
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        "<title>Gap Junction Thesis Showcase</title>",
        f"<style>{_GALLERY_CSS}</style></head><body>",
        "<div class='layout'>",
        "<nav><h2>Sections</h2>",
    ]
    for sid, label, _ in GALLERY_SECTIONS:
        if sid in sections:
            parts.append(f"<a href='#{sid}'>{label}</a>")
    parts.extend([
        "<div class='switch'>Story page<br>"
        "<a href='../gap_junction_path/index.html'>Guided story path →</a><br>"
        "<a href='../phylogenetic_story/index.html'>Phylogeny story (inx → panx → cnx) →</a><br>"
        "<a href='../family_comparison/index.html'>Innexin vs connexin →</a><br>"
        "<a href='../innexin_insights/index.html'>Innexin insights →</a><br>"
        "<a href='../connexin_insights/index.html'>Connexin insights →</a></div>",
        "</nav><main>",
        "<div class='hero'>",
        "<h1>Gap Junction Comparative Genomics</h1>",
        "<p class='subtitle'>Innexins · Connexins · Phylogeny · Domains · Synteny</p>",
        "<div class='stats'>",
        f"<span class='stat'><b>{stats['species']}</b> species</span>",
        f"<span class='stat'><b>{stats['candidates']}</b> candidates</span>",
        f"<span class='stat'><b>{stats['accepted']}</b> accepted</span>",
        f"<span class='stat'><b>{stats['curator_present_spp']}</b> curator present spp "
        f"(<b>≥{stats['curator_present_loci']}</b> loci, floor)</span>",
        "</div></div>",
    ])
    for sid, label, _intro in GALLERY_SECTIONS:
        if sid not in sections:
            continue
        parts.append(f"<section id='{sid}'><h2>{label}</h2><div class='grid'>")
        for fig_title, path in sections[sid]:
            rel = path.relative_to(out_dir).as_posix()
            parts.append(
                f"<div class='card'>"
                f"<img src='{rel}' alt='{_escape_html(fig_title)}' loading='lazy' "
                f"onclick=\"openModal(this.src)\">"
                f"<h3 class='peer-title'>{_escape_html(fig_title)}</h3></div>"
            )
        parts.append("</div></section>")
    parts.extend([
        "<div id='modal' class='modal' onclick='closeModal()'>",
        "<img id='modal-img' src='' alt='Enlarged figure'></div>",
        "<script>",
        "function openModal(src){document.getElementById('modal-img').src=src;",
        "document.getElementById('modal').classList.add('open')}",
        "function closeModal(){document.getElementById('modal').classList.remove('open')}",
        "document.addEventListener('keydown',e=>{if(e.key==='Escape')closeModal()})",
        "</script></main></div></body></html>",
    ])
    return "\n".join(parts)


def _render_presenter_gallery_html(
    sections: dict[str, list[tuple[str, str, Path]]],
    out_dir: Path,
) -> str:
    stats = _gallery_stats()
    lookup = _path_lookup(sections, out_dir)
    findings = _ordered_key_findings(sections, out_dir)
    parts: list[str] = [
        "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
        "<meta name='viewport' content='width=device-width,initial-scale=1'>",
        "<title>Key findings</title>",
        f"<style>{_GALLERY_CSS}</style></head><body>",
        "<div class='layout'>",
        "<nav><h2>Key findings</h2>",
    ]
    for i, finding in enumerate(findings, 1):
        rel = finding["rel"]
        headline = finding["headline"]
        if rel == "discovery/02_innexin_candidates.png":
            headline = (
                f"True rescues: {stats['curator_rescued_spp']} species "
                f"(≥{stats['curator_rescued_loci']} loci, floor); "
                f"+{stats['curator_firstpass_spp']} first-pass curator hits"
            )
        short = headline[:42] + ("…" if len(headline) > 42 else "")
        parts.append(
            f"<a href='#finding-{i}' title='{_escape_html(headline)}'>"
            f"{i}. {_escape_html(short)}</a>"
        )
    parts.append(
        "<div class='switch'>Full gallery<br>"
        "<a href='index.html'>Open showcase →</a></div>"
    )
    parts.extend([
        "</nav><main>",
        "<div class='hero presenter'>",
        "<h1>Key findings</h1>",
        "<p class='subtitle'>Figure notes · "
        "<a href='index.html' style='color:#fff'>index.html</a></p>",
        "<div class='stats'>",
        f"<span class='stat'><b>{len(findings)}</b> findings</span>",
        f"<span class='stat'><b>{stats['accepted']}</b> accepted candidates</span>",
        f"<span class='stat'><b>{stats['curator_rescued_spp']}</b> rescued spp · "
        f"<b>≥{stats['curator_rescued_loci']}</b> loci (floor)</span>",
        f"<span class='stat'><b>{stats['curator_firstpass_spp']}</b> first-pass spp · "
        f"<b>≥{stats['curator_firstpass_loci']}</b> loci</span>",
        "</div></div>",
        "<div class='guide'><h2>Order</h2><ol>",
    ])
    for i, finding in enumerate(findings, 1):
        rel = finding["rel"]
        headline = finding["headline"]
        if rel == "discovery/02_innexin_candidates.png":
            headline = (
                f"True rescues: {stats['curator_rescued_spp']} species "
                f"(≥{stats['curator_rescued_loci']} loci, floor); "
                f"+{stats['curator_firstpass_spp']} first-pass curator hits"
            )
        parts.append(
            f"<li><a href='#finding-{i}'>{_escape_html(headline)}</a></li>"
        )
    parts.extend([
        "</ol></div>",
        "<section id='findings'>",
        "<p class='section-intro'>Selected results in the same order as index.html.</p>",
        "<div class='findings-list'>",
    ])
    for i, finding in enumerate(findings, 1):
        rel = finding["rel"]
        fig_title, _ = lookup[rel]
        headline = finding["headline"]
        caption = finding.get("caption", PRESENTER_FIGURE_CAPTIONS.get(rel, ""))
        script = finding.get("script", PRESENTER_SPEAKER_SCRIPT.get(rel, ""))
        if rel == "discovery/02_innexin_candidates.png":
            headline = (
                f"True rescues: {stats['curator_rescued_spp']} species "
                f"(≥{stats['curator_rescued_loci']} loci, floor); "
                f"+{stats['curator_firstpass_spp']} first-pass curator hits"
            )
            script = (
                f"Correct framing: {stats['curator_rescued_spp']} true rescues "
                f"(prior discovery accepted=0; ≥{stats['curator_rescued_loci']} loci). "
                f"Separately, {stats['curator_firstpass_spp']} species were never run through "
                f"discovery (not_run) and first show curator hits "
                f"(≥{stats['curator_firstpass_loci']} loci). "
                f"Do not call the whole set 'rescued'."
            )
        cap_html = f"<p class='caption'>{_escape_html(caption)}</p>" if caption else ""
        script_html = f"<p class='script'>{_escape_html(script)}</p>" if script else ""
        parts.append(
            f"<article class='finding-card key' id='finding-{i}'>"
            f"<span class='badge'>Finding {i}</span>"
            f"<p class='finding-headline'>{_escape_html(headline)}</p>"
            f"<p class='caption' style='margin-top:0'><em>{_escape_html(fig_title)}</em></p>"
            f"<img src='{rel}' alt='{_escape_html(fig_title)}' loading='lazy' "
            f"onclick=\"openModal(this.src)\">"
            f"{cap_html}{script_html}</article>"
        )
    parts.extend([
        "</div></section>",
        "<div id='modal' class='modal' onclick='closeModal()'>",
        "<img id='modal-img' src='' alt='Enlarged figure'></div>",
        "<script>",
        "function openModal(src){document.getElementById('modal-img').src=src;",
        "document.getElementById('modal').classList.add('open')}",
        "function closeModal(){document.getElementById('modal').classList.remove('open')}",
        "document.addEventListener('keydown',e=>{if(e.key==='Escape')closeModal()})",
        "</script></main></div></body></html>",
    ])
    return "\n".join(parts)


def write_gallery_index(plot_records: list[tuple[str, str, Path]], out: Path) -> None:
    out_dir = out.parent
    sections = _prepare_gallery_sections(plot_records, out_dir)
    out.write_text(_render_peer_gallery_html(sections, out_dir), encoding="utf-8")
    presenter = out_dir / "presenter.html"
    presenter.write_text(_render_presenter_gallery_html(sections, out_dir), encoding="utf-8")


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def build_showcase(out_dir: Path) -> list[tuple[str, str, Path]]:
    tables = load_tables()
    ref_inx = load_reference_lengths("innexin")
    ref_cnx = load_reference_lengths("connexin")
    ref_all = pd.concat([ref_inx, ref_cnx], ignore_index=True)
    records: list[tuple[str, str, Path]] = []

    def add(section: str, title: str, rel: str, fn) -> None:
        path = out_dir / rel
        fn(path)
        records.append((section, title, path))

    # Overview
    add("overview", "Species panel composition", "overview/01_panel_composition.png",
        lambda p: plot_species_panel_composition(tables["panel"], p))
    add("overview", "Discovery funnel", "overview/02_discovery_funnel.png",
        lambda p: plot_discovery_funnel(tables["candidates"], p))
    add("overview", "Annotation status", "overview/04_annotation_status.png",
        lambda p: plot_annotation_status(tables["thesis"], p))

    # Innexin
    add("innexin", "Innexin protein length by clade", "innexin/01_length_by_clade.png",
        lambda p: plot_clade_protein_lengths(ref_inx, "innexin", p))
    add("innexin", "Innexin exon count by organism", "innexin/02_exon_by_organism.png",
        lambda p: plot_exon_count_by_organism(tables["gene_summary"], "innexin", p))
    add("innexin", "Innexin gene span vs exons", "innexin/03_gene_span_vs_exons.png",
        lambda p: plot_gene_span_vs_exons(tables["gene_summary"], "innexin", p))

    # Connexin
    add("connexin", "Connexin protein length by clade", "connexin/01_length_by_clade.png",
        lambda p: plot_clade_protein_lengths(ref_cnx, "connexin", p))
    add("connexin", "Connexin exon count by organism", "connexin/02_exon_by_organism.png",
        lambda p: plot_exon_count_by_organism(tables["gene_summary"], "connexin", p))
    add("connexin", "Connexin gene span vs exons", "connexin/03_gene_span_vs_exons.png",
        lambda p: plot_gene_span_vs_exons(tables["gene_summary"], "connexin", p))

    # Comparison
    add("comparison", "Protein length innexin vs connexin", "comparison/01_protein_length_boxplot.png",
        lambda p: plot_protein_length_by_family(ref_all, p))
    add("comparison", "Discovery scatter length vs TM", "comparison/02_discovery_scatter.png",
        lambda p: plot_discovery_scatter(tables["candidates"], p))
    add("comparison", "Rank score distribution", "comparison/03_rank_scores.png",
        lambda p: plot_rank_score_distribution(tables["candidates"], p))
    add("comparison", "Exon count comparison", "comparison/04_exon_count_boxplot.png",
        lambda p: _exon_family_boxplot(tables["gene_summary"], p))
    add("comparison", "Exon vs TM in discovery", "comparison/05_exon_vs_tm.png",
        lambda p: plot_exon_tm_comparison(tables["candidates"], p))

    add("comparison", "Acceptance rates", "comparison/06_acceptance_rates.png",
        lambda p: plot_acceptance_rates(tables["candidates"], p))
    add("comparison", "Protein length histograms", "comparison/07_length_histograms.png",
        lambda p: plot_protein_length_histogram(ref_all, p))
    add("comparison", "Intron architecture", "comparison/08_intron_architecture.png",
        lambda p: plot_intron_architecture(tables["transcript_summary"], p))

    # Innexin extras
    add("innexin", "Dmel innexin paralog lengths", "innexin/04_dmel_paralog_lengths.png",
        lambda p: plot_dmel_innexin_lengths(ref_inx, p))
    add("innexin", "Innexin genes per species", "innexin/05_genes_per_species.png",
        lambda p: plot_genes_per_species(tables["gene_summary"], p))

    # Connexin extras
    add("connexin", "Human connexin landscape", "connexin/04_human_landscape.png",
        lambda p: plot_human_connexin_landscape(ref_cnx, p))
    add("connexin", "Connexin genes per species", "connexin/05_genes_per_species.png",
        lambda p: plot_genes_per_species(tables["gene_summary"], p))

    # Domains
    domain_dir = out_dir / "domains"
    plot_domain_coverage_maps(tables["features"], tables["structure"], domain_dir)
    for name, title in [
        ("domain_mv_pfam_connexin.png", "Connexin Pfam monument valley"),
        ("domain_mv_tm_topology.png", "TM topology coverage"),
        ("domain_mv_disorder.png", "Disorder coverage"),
        ("domain_mv_coiled_coil.png", "Coiled-coil coverage"),
    ]:
        p = domain_dir / name
        if p.exists():
            records.append(("domains", title, p))
    add("domains", "Feature fraction heatmap", "domains/05_feature_fraction_heatmap.png",
        lambda p: plot_feature_fraction_bars(tables["structure"], p))

    # Discovery
    add("discovery", "Validation by species", "discovery/01_validation_by_species.png",
        lambda p: plot_validation_by_species(tables["validation"], p))
    add("discovery", "Innexin-only candidates", "discovery/02_innexin_candidates.png",
        lambda p: _innexin_candidate_panel(tables["candidates"], p))

    # Phylogeny
    tree = PHYLO_DIR / "gap_junction_phylogeny_input.aln.treefile"
    add("phylogeny", "Gap junction phylogeny (family-colored)", "phylogeny/01_tree_colored.png",
        lambda p: plot_phylogeny_colored(tree, p, "Gap junction phylogeny — references + discovery"))
    add("phylogeny", "Phylogenetic distance heatmap", "phylogeny/02_distance_heatmap.png",
        lambda p: plot_phylo_distance_heatmap(tree, p))
    inx_tree = PHYLO_DIR / "innexin_reference_panel.treefile"
    add("phylogeny", "Innexin reference panel tree", "phylogeny/03_innexin_reference_tree.png",
        lambda p: plot_phylogeny_colored(inx_tree, p, "Innexin reference panel phylogeny"))

    # Synteny
    add("synteny", "SynVoy job summary", "synteny/01_job_summary.png",
        lambda p: plot_synteny_job_summary(tables["synteny_jobs"], p))
    add("synteny", "SynVoy annotation yield", "synteny/02_annotation_yield.png",
        lambda p: plot_synteny_annotations(tables["synteny_jobs"], p))
    add("synteny", "Reference panel SynVoy status", "synteny/03_refpanel_status.png",
        lambda p: _refpanel_synteny_status(p))
    add("synteny", "Case 1 Dmel cluster sharing", "synteny/04_case1_cluster_sharing.png",
        lambda p: plot_case1_synteny_sharing(p))
    add("synteny", "GOI annotation heatmap", "synteny/05_goi_heatmap.png",
        lambda p: plot_synteny_goi_heatmap(p))
    add("synteny", "Flanking vs GOI counts", "synteny/06_flanking_goi_scatter.png",
        lambda p: plot_flanking_gene_counts(p))

    add("cross_clade", "Case 3 insect vs nematode", "cross_clade/01_case3_comparison.png",
        lambda p: plot_case3_cross_clade(p))
    add("cross_clade", "Inx2 insect vs nematode bars", "cross_clade/02_inx2_architecture.png",
        lambda p: plot_inx2_vs_inx2_nematode_insect(tables["gene_summary"], ref_inx, p))

    add("ortholog", "Case 2 ortholog support", "ortholog/01_case2_support.png",
        lambda p: plot_case2_ortholog_support(p))
    add("ortholog", "Inx2/shakB ortholog panel", "ortholog/02_inx2_panel.png",
        lambda p: plot_inx2_ortholog_panel(tables["gene_summary"], ref_inx, p))
    add("ortholog", "shakB insect comparison", "ortholog/03_shakB_insects.png",
        lambda p: plot_shakB_insect_comparison(ref_inx, p))

    add("connexin", "GJA1 ortholog trajectory", "connexin/06_gja1_trajectory.png",
        lambda p: plot_gja1_ortholog_trajectory(ref_cnx, p))

    add("comparison", "Reference vs discovery lengths", "comparison/09_ref_vs_discovery.png",
        lambda p: plot_reference_vs_discovery(tables["candidates"], ref_all, p))
    add("comparison", "Rank category breakdown", "comparison/10_rank_breakdown.png",
        lambda p: plot_rank_category_breakdown(tables["candidates"], p))
    add("comparison", "Exon length distribution", "comparison/11_exon_length_dist.png",
        lambda p: plot_exon_length_distribution(p))
    add("comparison", "Transcript complexity", "comparison/12_transcript_complexity.png",
        lambda p: plot_transcript_complexity(tables["gene_summary"], p))

    add("innexin", "Nematode innexin panel", "innexin/06_nematode_panel.png",
        lambda p: plot_nematode_innexin_panel(ref_inx, p))

    add("domains", "TM helices vs length", "domains/06_tm_vs_length.png",
        lambda p: plot_tm_helix_reference(tables["structure"], ref_all, p))
    add("domains", "Secondary structure fractions", "domains/07_secstruct_fraction.png",
        lambda p: plot_secondary_structure_fraction(tables["structure"], ref_all, p))

    add("discovery", "Connexin batch comparison", "discovery/03_connexin_batches.png",
        lambda p: plot_connexin_batch_discovery(tables["validation"], p))

    write_gallery_index(records, out_dir / "index.html")
    manifest = out_dir / "plot_manifest.csv"
    with manifest.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(["section", "title", "path"])
        for section, title, path in records:
            if path.exists():
                writer.writerow([section, title, str(path)])
    return records


def _exon_family_boxplot(gene_summary: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    data = [gene_summary[gene_summary["family"] == f]["mean_exon_count"].values for f in ("innexin", "connexin")]
    bp = ax.boxplot(data, tick_labels=["Innexin", "Connexin"], patch_artist=True)
    for patch, fam in zip(bp["boxes"], ("innexin", "connexin")):
        patch.set_facecolor(FAMILY_COLORS[fam])
        patch.set_alpha(0.6)
    ax.set_ylabel("Mean exon count per gene")
    ax.set_title("Exon architecture: innexin vs connexin")
    _save(fig, out)


def _innexin_candidate_panel(candidates: pd.DataFrame, out: Path) -> None:
    sub = candidates[
        (candidates["family"] == "innexin") & (candidates["rank_category"] != "rejected_false_positive")
    ].copy()
    if sub.empty:
        fig, ax = plt.subplots(figsize=(8, 3))
        ax.text(0.5, 0.5, "No accepted innexin candidates", ha="center", va="center")
        ax.axis("off")
        _save(fig, out)
        return
    # Aggregate by organism when many curator loci (readable overview)
    if len(sub) > 25:
        grouped = (
            sub.groupby("organism", as_index=False)
            .agg(
                n_loci=("candidate_id", "count"),
                best_aa=("protein_length", "max"),
                curator=("discovery_batch", lambda s: int((s == "gene_curator_probe").any())),
            )
            .sort_values("best_aa")
        )
        fig, ax = plt.subplots(figsize=(10, max(5, 0.35 * len(grouped))))
        colors = ["#27AE60" if c else "#F39C12" for c in grouped["curator"]]
        ax.barh(grouped["organism"], grouped["best_aa"], color=colors, alpha=0.85)
        for y, (_, row) in enumerate(grouped.iterrows()):
            ax.text(row["best_aa"] + 5, y, f"n={int(row['n_loci'])}", va="center", fontsize=8)
        ax.set_xlabel("Longest accepted product (aa)")
        ax.set_title(
            f"Accepted innexin candidates by species (n={len(sub)} loci, {len(grouped)} species)"
        )
    else:
        fig, ax = plt.subplots(figsize=(9, max(4, 0.35 * len(sub))))
        colors = [
            "#27AE60" if ("high" in r or "curator_present" in r) else "#F39C12"
            for r in sub["rank_category"]
        ]
        ax.barh(
            sub["organism"] + " | " + sub["candidate_id"],
            sub["protein_length"],
            color=colors,
            alpha=0.85,
        )
        ax.set_xlabel("Protein length (aa)")
        ax.set_title("Accepted innexin discovery candidates")
    _save(fig, out)


def plot_protein_length_histogram(ref_all: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.5), sharey=True)
    for ax, family in zip(axes, ("innexin", "connexin")):
        sub = ref_all[ref_all["family"] == family]
        ax.hist(sub["protein_length"], bins=20, color=FAMILY_COLORS[family], alpha=0.75, edgecolor="white")
        ax.set_xlabel("Protein length (aa)")
        ax.set_title(f"{family.capitalize()} reference proteins (n={len(sub)})")
    axes[0].set_ylabel("Count")
    _save(fig, out)


def plot_genes_per_species(gene_summary: pd.DataFrame, out: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    for ax, family in zip(axes, ("innexin", "connexin")):
        sub = gene_summary[gene_summary["family"] == family].groupby("organism").size().sort_values()
        ax.barh(sub.index, sub.values, color=FAMILY_COLORS[family], alpha=0.85)
        ax.set_xlabel("Reference genes in panel")
        ax.set_title(f"{family.capitalize()} genes per species")
    _save(fig, out)


def plot_intron_architecture(transcript_summary: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(8, 5))
    for family in ("innexin", "connexin"):
        sub = transcript_summary[transcript_summary["family"] == family]
        ax.scatter(
            sub["exon_count"],
            sub["mean_intron_length"],
            c=FAMILY_COLORS[family],
            alpha=0.35,
            s=25,
            label=family,
        )
    ax.set_xlabel("Exon count")
    ax.set_ylabel("Mean intron length (bp)")
    ax.set_yscale("log")
    ax.set_title("Intron architecture: innexin vs connexin")
    ax.legend()
    _save(fig, out)


def plot_case1_synteny_sharing(out: Path) -> None:
    """Two genomic clusters on Dmel X — not one bar per paralog×region."""
    csv_path = RESULTS_DIR / "synteny_case_studies" / "case1_dmel_innexin_cluster" / "case1_region_summary.csv"
    # Prefer curated two-cluster summary if present
    curated = RESULTS_DIR / "synteny_case_studies" / "case1_dmel_innexin_cluster" / "case1_two_clusters.csv"
    if curated.exists():
        df = pd.read_csv(curated)
        fig, ax = plt.subplots(figsize=(9, 3.8))
        y = range(len(df))
        ax.barh(y, df["shared_flanking"], color="#2E86AB", alpha=0.9, label="Shared flanking genes")
        ax.barh(
            y,
            df["unique_flanking"],
            left=df["shared_flanking"],
            color="#E67E22",
            alpha=0.85,
            label="Unique flanking",
        )
        labels = [
            f"{row.cluster_label}\n{row.innexin_genes}"
            for row in df.itertuples()
        ]
        ax.set_yticks(list(y))
        ax.set_yticklabels(labels, fontsize=9)
        ax.set_xlabel("Flanking gene count")
        ax.set_title("D. melanogaster: two innexin genomic clusters (X chromosome)")
        ax.legend(fontsize=8, loc="lower right")
        _save(fig, out)
        return

    if not csv_path.exists():
        return
    df = pd.read_csv(csv_path)
    # Aggregate SynVoy case1 rows into the two genomic regions
    agg = (
        df.groupby("region", as_index=False)
        .agg(
            shared_flanking=("shared_with_all_paralogs", "max"),
            unique_flanking=("unique_to_query", "sum"),
            genes=("query_gene", lambda s: " · ".join(sorted(set(s)))),
        )
    )
    region_map = {
        "proximal_cluster_6.8-7.1Mb": ("Proximal X ~6.9–7.0 Mb", "ogre · Inx7 · Inx2"),
        "distal_cluster_20.5-21.0Mb": ("Distal X ~20.7–20.9 Mb", "shakB (+ SynVoy GOI models)"),
    }
    fig, ax = plt.subplots(figsize=(9, 3.8))
    y = range(len(agg))
    ax.barh(y, agg["shared_flanking"], color="#2E86AB", alpha=0.9, label="Shared flanking genes")
    ax.barh(
        y,
        agg["unique_flanking"],
        left=agg["shared_flanking"],
        color="#E67E22",
        alpha=0.85,
        label="Unique flanking",
    )
    labels = []
    for row in agg.itertuples():
        title, genes = region_map.get(row.region, (row.region, row.genes))
        labels.append(f"{title}\n{genes}")
    ax.set_yticks(list(y))
    ax.set_yticklabels(labels, fontsize=9)
    ax.set_xlabel("Flanking gene count")
    ax.set_title("D. melanogaster: two innexin genomic clusters (not five separate loci)")
    ax.legend(fontsize=8, loc="lower right")
    _save(fig, out)


def plot_dmel_innexin_lengths(ref_inx: pd.DataFrame, out: Path) -> None:
    sub = ref_inx[ref_inx["organism"] == "Drosophila melanogaster"].sort_values("protein_length")
    fig, ax = plt.subplots(figsize=(8, 5))
    ax.barh(sub["gene_symbol"], sub["protein_length"], color="#2E86AB", alpha=0.85)
    ax.set_xlabel("Protein length (aa)")
    ax.set_title("D. melanogaster innexin paralog lengths")
    _save(fig, out)


def plot_human_connexin_landscape(ref_cnx: pd.DataFrame, out: Path) -> None:
    sub = ref_cnx[ref_cnx["organism"] == "Homo sapiens"].sort_values("protein_length")
    fig, ax = plt.subplots(figsize=(9, 6))
    ax.barh(sub["gene_symbol"], sub["protein_length"], color="#E67E22", alpha=0.85)
    ax.set_xlabel("Protein length (aa)")
    ax.set_title("Human connexin gene family — protein lengths")
    _save(fig, out)


def plot_acceptance_rates(candidates: pd.DataFrame, out: Path) -> None:
    fig, ax = plt.subplots(figsize=(7, 5))
    rows = []
    for family in ("innexin", "connexin"):
        sub = candidates[candidates["family"] == family]
        rows.append({
            "family": family,
            "accepted": len(sub[sub["rank_category"] != "rejected_false_positive"]),
            "rejected": len(sub[sub["rank_category"] == "rejected_false_positive"]),
        })
    df = pd.DataFrame(rows)
    x = np.arange(2)
    ax.bar(x, df["accepted"], color="#27AE60", label="Accepted")
    ax.bar(x, df["rejected"], bottom=df["accepted"], color="#BDC3C7", label="Rejected")
    ax.set_xticks(x)
    ax.set_xticklabels(df["family"])
    ax.set_ylabel("Candidates")
    ax.set_title("Discovery acceptance rate by family")
    ax.legend()
    _save(fig, out)


def _refpanel_synteny_status(out: Path) -> None:
    jobs = [
        ("Inx2", "refpanel_dmel_Inx2"),
        ("shakB", "refpanel_dmel_shakB"),
        ("inx-2", "refpanel_cele_inx-2"),
    ]
    rows = []
    for gene, dirname in jobs:
        report = SYNVOY_RESULTS / dirname / "synvoy_report.json"
        plots = len(list((SYNVOY_RESULTS / dirname).glob("synteny_block_*_synteny_plot.html"))) if (SYNVOY_RESULTS / dirname).exists() else 0
        status = "complete" if report.exists() else "missing"
        rows.append((gene, plots, status))
    fig, ax = plt.subplots(figsize=(7, 4))
    genes, plot_counts, _ = zip(*rows)
    colors = ["#27AE60" if s == "complete" else "#E74C3C" for _, _, s in rows]
    ax.bar(genes, plot_counts, color=colors, alpha=0.85)
    ax.set_ylabel("Synteny plots")
    ax.set_title("Innexin reference panel SynVoy (completed)")
    _save(fig, out)


def main() -> None:
    parser = argparse.ArgumentParser(description="Build scientific showcase plots")
    parser.add_argument("--outdir", type=Path, default=SHOWCASE_DIR)
    parser.add_argument(
        "--html-only",
        action="store_true",
        help="Rebuild index.html / presenter.html from existing PNGs (no plot regen)",
    )
    parser.add_argument(
        "--refresh-discovery",
        action="store_true",
        help="Also regenerate discovery/overview plots that depend on candidate metrics",
    )
    args = parser.parse_args()
    out_dir = args.outdir
    out_dir.mkdir(parents=True, exist_ok=True)

    if args.html_only or args.refresh_discovery:
        tables = load_tables()
        if args.refresh_discovery:
            (out_dir / "overview").mkdir(parents=True, exist_ok=True)
            (out_dir / "discovery").mkdir(parents=True, exist_ok=True)
            (out_dir / "comparison").mkdir(parents=True, exist_ok=True)
            plot_species_panel_composition(tables["panel"], out_dir / "overview/01_panel_composition.png")
            plot_discovery_funnel(tables["candidates"], out_dir / "overview/02_discovery_funnel.png")
            plot_discovery_scatter(tables["candidates"], out_dir / "comparison/02_discovery_scatter.png")
            plot_rank_score_distribution(tables["candidates"], out_dir / "comparison/03_rank_scores.png")
            plot_acceptance_rates(tables["candidates"], out_dir / "comparison/06_acceptance_rates.png")
            plot_rank_category_breakdown(tables["candidates"], out_dir / "comparison/10_rank_breakdown.png")
            _innexin_candidate_panel(tables["candidates"], out_dir / "discovery/02_innexin_candidates.png")

        records = _collect_existing_plot_records(out_dir)
        write_gallery_index(records, out_dir / "index.html")
        print(f"HTML refreshed from {len(records)} existing plots")
        print(f"Gallery:         {out_dir / 'index.html'}")
        print(f"Presenter page:  {out_dir / 'presenter.html'}")
        return

    records = build_showcase(out_dir)
    n = sum(1 for _, _, p in records if p.exists())
    print(f"Showcase complete: {n} plots in {out_dir}")
    print(f"Gallery:         {out_dir / 'index.html'}")
    print(f"Presenter page:  {out_dir / 'presenter.html'}")


def _collect_existing_plot_records(out_dir: Path) -> list[tuple[str, str, Path]]:
    """Walk showcase section folders and rebuild (section, title, path) records."""
    section_map = {sid: label for sid, label, _ in GALLERY_SECTIONS}
    records: list[tuple[str, str, Path]] = []
    for sid in section_map:
        folder = out_dir / sid
        if not folder.is_dir():
            continue
        for path in sorted(folder.rglob("*.png")):
            title = _clean_plot_title(path.stem)
            records.append((sid, title, path))
    # imported nests
    for path in sorted((out_dir / "imported").rglob("*.png")) if (out_dir / "imported").is_dir() else []:
        parent = path.parent.name
        sid = "exon_maps" if parent == "exon_maps" else "protein_features" if parent == "protein_features" else "domains"
        records.append((sid, _clean_plot_title(path.stem), path))
    return records


if __name__ == "__main__":
    main()
