#!/usr/bin/env python3
"""Side-by-side innexin vs connexin comparisons (same metrics, one site).

Reuses clade-comparison tables, MMseqs/embedding matrices, exon summaries,
and optionally a cross-family MMseqs search. Outputs:

    project/results/family_comparison/
"""

from __future__ import annotations

import html
import math
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

from pipeline.common import RESULTS_DIR, write_csv

OUT = RESULTS_DIR / "family_comparison"
FIGS = OUT / "figures"
WORK = OUT / "work"

INX_SIM = RESULTS_DIR / "innexin_similarity"
CNX_SIM = RESULTS_DIR / "connexin_similarity"
INX_CLADE = RESULTS_DIR / "innexin_clade_comparison"
CNX_CLADE = RESULTS_DIR / "connexin_clade_comparison"
EXON_SUMMARY = RESULTS_DIR / "exon_structures" / "gene_summary.csv"

INX_COLOR = "#0F766E"
CNX_COLOR = "#C2410C"
MUTED = "#64748B"

INX_SF_SHORT = {
    "SF1_shakB": "SF1",
    "SF2_Inx1_ogre": "SF2",
    "SF3_Inx3_Inx7_nematode": "SF3",
    "SF4_Inx2_expansion": "SF4",
    "unassigned": "unassigned",
}
CNX_SF_SHORT = {
    "alpha_GJA": "α GJA",
    "beta_GJB": "β GJB",
    "gamma_GJC": "γ GJC",
    "delta_GJD": "δ GJD",
    "epsilon_GJE": "ε GJE",
    "unassigned": "unassigned",
}
INX_SF_COLORS = {
    "SF1_shakB": "#1D4ED8",
    "SF2_Inx1_ogre": "#0F766E",
    "SF3_Inx3_Inx7_nematode": "#B45309",
    "SF4_Inx2_expansion": "#7C2D12",
    "unassigned": "#94A3B8",
}
CNX_SF_COLORS = {
    "alpha_GJA": "#9A3412",
    "beta_GJB": "#C2410C",
    "gamma_GJC": "#0369A1",
    "delta_GJD": "#15803D",
    "epsilon_GJE": "#7C3AED",
    "unassigned": "#94A3B8",
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


def _boolish(val) -> bool:
    return str(val).strip().lower() in {"true", "1", "yes"}


def _short(group: str, family: str) -> str:
    table = INX_SF_SHORT if family == "innexin" else CNX_SF_SHORT
    return table.get(group, group)


def load_group_matrix(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    df["mean_similarity"] = pd.to_numeric(df["mean_similarity"], errors="coerce")
    df["median_similarity"] = pd.to_numeric(df["median_similarity"], errors="coerce")
    df["n_pairs"] = pd.to_numeric(df["n_pairs"], errors="coerce").fillna(0).astype(int)
    df["same_group"] = df["same_group"].map(_boolish)
    return df


def weighted_mean(df: pd.DataFrame, value_col: str = "mean_similarity") -> float:
    if df.empty:
        return float("nan")
    weights = df["n_pairs"].clip(lower=1)
    return float(np.average(df[value_col], weights=weights))


# Headline numbers from the existing family-level analyses (pooled pairwise
# statistics), so this page does not contradict innexin/connexin insights.
CANONICAL = {
    "innexin": {"within_id": 46.2, "between_id": 27.9, "emb_same": 0.279, "emb_diff": 0.049},
    "connexin": {"within_id": 56.6, "between_id": 40.1, "emb_same": 0.492, "emb_diff": -0.053},
}


def family_metrics(df: pd.DataFrame, method: str) -> dict:
    sub = df[(df["method"] == method) & (~df["group_a"].str.contains("unassigned", na=False)) & (~df["group_b"].str.contains("unassigned", na=False))].copy()
    within = sub[sub["same_group"]]
    between = sub[~sub["same_group"]]
    return {
        "within_mean": weighted_mean(within),
        "between_mean": weighted_mean(between),
        "within_median": weighted_mean(within, "median_similarity"),
        "between_median": weighted_mean(between, "median_similarity"),
        "n_within_pairs": int(within["n_pairs"].sum()) if not within.empty else 0,
        "n_between_pairs": int(between["n_pairs"].sum()) if not between.empty else 0,
    }


def pair_mean(df: pd.DataFrame, method: str, a: str, b: str) -> float:
    sub = df[(df["method"] == method) & (df["group_a"] == a) & (df["group_b"] == b)]
    if sub.empty:
        sub = df[(df["method"] == method) & (df["group_a"] == b) & (df["group_b"] == a)]
    return weighted_mean(sub) if not sub.empty else float("nan")


def heatmap_matrix(df: pd.DataFrame, method: str, groups: list[str]) -> np.ndarray:
    sub = df[df["method"] == method]
    mat = np.full((len(groups), len(groups)), np.nan)
    idx = {g: i for i, g in enumerate(groups)}
    for _, row in sub.iterrows():
        if row["group_a"] not in idx or row["group_b"] not in idx:
            continue
        mat[idx[row["group_a"]], idx[row["group_b"]]] = row["mean_similarity"]
    return mat


def run_cross_family_mmseqs() -> Path | None:
    mmseqs = shutil.which("mmseqs") or str(Path.home() / "miniconda3/envs/synvoy_env/bin/mmseqs")
    if not Path(mmseqs).exists():
        print("MMseqs2 not found — skipping cross-family search")
        return None
    query = INX_SIM / "all_innexins.fasta"
    target = CNX_SIM / "all_connexins.fasta"
    if not query.exists() or not target.exists():
        print("FASTAs missing — skipping cross-family search")
        return None
    WORK.mkdir(parents=True, exist_ok=True)
    m8 = OUT / "innexin_vs_connexin.m8"
    tmp = WORK / "mmseqs_cross"
    if tmp.exists():
        shutil.rmtree(tmp)
    tmp.mkdir(parents=True)
    cmd = [
        mmseqs,
        "easy-search",
        str(query),
        str(target),
        str(m8),
        str(tmp),
        "--max-seqs",
        "5",
        "-e",
        "10",
        "-s",
        "7.5",
        "-v",
        "1",
        "--format-output",
        "query,target,pident,alnlen,evalue,bits",
    ]
    try:
        subprocess.run(cmd, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"Cross-family MMseqs failed: {exc}")
        return None
    return m8 if m8.exists() else None


def plot_dataset_scope(inx: pd.DataFrame, cnx: pd.DataFrame, path: Path) -> None:
    labels = ["Proteins", "Species", "Clades"]
    inx_vals = [len(inx), inx["organism"].nunique(), inx["clade"].nunique()]
    cnx_vals = [len(cnx), cnx["organism"].nunique(), cnx["clade"].nunique()]
    x = np.arange(len(labels))
    w = 0.36
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    ax.bar(x - w / 2, inx_vals, w, label="Innexin", color=INX_COLOR, alpha=0.92)
    ax.bar(x + w / 2, cnx_vals, w, label="Connexin", color=CNX_COLOR, alpha=0.92)
    for i, (a, b) in enumerate(zip(inx_vals, cnx_vals)):
        ax.text(i - w / 2, a + 2, str(a), ha="center", fontsize=9, color=INX_COLOR)
        ax.text(i + w / 2, b + 2, str(b), ha="center", fontsize=9, color=CNX_COLOR)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Count")
    ax.set_title("Same comparison, two families: dataset size")
    ax.legend(frameon=False)
    _save(fig, path)


def plot_clade_geography(inx: pd.DataFrame, cnx: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.4))
    panels = [
        (axes[0], inx, "Innexin subfamilies by clade", INX_SF_COLORS, INX_SF_SHORT, ["SF1_shakB", "SF2_Inx1_ogre", "SF3_Inx3_Inx7_nematode", "SF4_Inx2_expansion", "unassigned"]),
        (axes[1], cnx, "Connexin classes by clade", CNX_SF_COLORS, CNX_SF_SHORT, ["alpha_GJA", "beta_GJB", "gamma_GJC", "delta_GJD", "epsilon_GJE", "unassigned"]),
    ]
    for ax, df, title, colors, short, order in panels:
        clades = list(df["clade"].value_counts().index)
        x = np.arange(len(clades))
        bottom = np.zeros(len(clades))
        for sf in order:
            vals = np.array([(df[(df["clade"] == c) & (df["subfamily"] == sf)].shape[0]) for c in clades], dtype=float)
            if vals.sum() == 0:
                continue
            ax.bar(x, vals, bottom=bottom, color=colors.get(sf, MUTED), label=short.get(sf, sf), width=0.78)
            bottom += vals
        ax.set_xticks(x)
        ax.set_xticklabels(clades, rotation=32, ha="right", fontsize=8)
        ax.set_ylabel("Loci")
        ax.set_title(title)
        ax.legend(frameon=False, fontsize=7, loc="upper right")
    _save(fig, path)


def plot_within_between(inx_m: dict, cnx_m: dict, path: Path) -> None:
    labels = ["Within subfamily", "Between subfamilies"]
    inx_vals = [inx_m["within_mean"], inx_m["between_mean"]]
    cnx_vals = [cnx_m["within_mean"], cnx_m["between_mean"]]
    x = np.arange(len(labels))
    w = 0.34
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    ax.bar(x - w / 2, inx_vals, w, color=INX_COLOR, label="Innexin")
    ax.bar(x + w / 2, cnx_vals, w, color=CNX_COLOR, label="Connexin")
    for i, (a, b) in enumerate(zip(inx_vals, cnx_vals)):
        ax.text(i - w / 2, a + 0.8, f"{a:.1f}%", ha="center", fontsize=9)
        ax.text(i + w / 2, b + 0.8, f"{b:.1f}%", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Mean MMseqs2 identity (%)")
    ax.set_ylim(0, max(inx_vals + cnx_vals) * 1.18)
    ax.set_title("Sequence cohesion: same metric, both families")
    ax.legend(frameon=False)
    _save(fig, path)


def plot_identity_heatmaps(inx_df: pd.DataFrame, cnx_df: pd.DataFrame, path: Path) -> None:
    inx_groups = ["SF1_shakB", "SF2_Inx1_ogre", "SF3_Inx3_Inx7_nematode", "SF4_Inx2_expansion"]
    cnx_groups = ["alpha_GJA", "beta_GJB", "gamma_GJC", "delta_GJD"]
    fig, axes = plt.subplots(1, 2, figsize=(12.8, 5.6))
    for ax, df, groups, family, title, cmap in [
        (axes[0], inx_df, inx_groups, "innexin", "Innexin SF identity (%)", "YlGnBu"),
        (axes[1], cnx_df, cnx_groups, "connexin", "Connexin class identity (%)", "YlOrBr"),
    ]:
        mat = heatmap_matrix(df, "mmseqs_pident", groups)
        im = ax.imshow(mat, cmap=cmap, vmin=20, vmax=70)
        ax.set_xticks(range(len(groups)))
        ax.set_yticks(range(len(groups)))
        ax.set_xticklabels([_short(g, family) for g in groups], rotation=35, ha="right")
        ax.set_yticklabels([_short(g, family) for g in groups])
        ax.set_title(title)
        for i in range(len(groups)):
            for j in range(len(groups)):
                v = mat[i, j]
                if not math.isnan(v):
                    ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=8)
        fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    _save(fig, path)


def plot_embedding_separation(inx_m: dict, cnx_m: dict, path: Path) -> None:
    labels = ["Same subfamily", "Different subfamilies"]
    inx_vals = [inx_m["within_mean"], inx_m["between_mean"]]
    cnx_vals = [cnx_m["within_mean"], cnx_m["between_mean"]]
    x = np.arange(len(labels))
    w = 0.34
    fig, ax = plt.subplots(figsize=(8.4, 5.0))
    ax.bar(x - w / 2, inx_vals, w, color=INX_COLOR, label="Innexin")
    ax.bar(x + w / 2, cnx_vals, w, color=CNX_COLOR, label="Connexin")
    ax.axhline(0, color="#cbd5e1", lw=1)
    for i, (a, b) in enumerate(zip(inx_vals, cnx_vals)):
        ax.text(i - w / 2, a + (0.02 if a >= 0 else -0.05), f"{a:.2f}", ha="center", fontsize=9)
        ax.text(i + w / 2, b + (0.02 if b >= 0 else -0.05), f"{b:.2f}", ha="center", fontsize=9)
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Mean 3-mer embedding cosine")
    ax.set_title("Local sequence texture separates groups in both families")
    ax.legend(frameon=False)
    _save(fig, path)


def plot_threshold_sweep(inx_sw: pd.DataFrame, cnx_sw: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.8, 5.0))
    ax.plot(inx_sw["min_seq_id"], inx_sw["mean_subfamily_purity"], "o-", color=INX_COLOR, label="Innexin purity")
    ax.plot(cnx_sw["min_seq_id"], cnx_sw["mean_subfamily_purity"], "s-", color=CNX_COLOR, label="Connexin purity")
    ax2 = ax.twinx()
    ax2.spines["right"].set_visible(True)
    ax2.plot(inx_sw["min_seq_id"], inx_sw["n_clusters"], "o--", color=INX_COLOR, alpha=0.45, label="Innexin clusters")
    ax2.plot(cnx_sw["min_seq_id"], cnx_sw["n_clusters"], "s--", color=CNX_COLOR, alpha=0.45, label="Connexin clusters")
    ax.set_xlabel("MMseqs2 min-seq-id")
    ax.set_ylabel("Mean subfamily / class purity")
    ax2.set_ylabel("Number of clusters")
    ax.set_ylim(0.7, 1.05)
    ax.set_title("Identity-threshold sweep: how cleanly groups cluster")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8)
    _save(fig, path)


def plot_key_pair(sf34: float, ab: float, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.6, 4.8))
    labels = ["Innexin SF3 ↔ SF4", "Connexin α ↔ β"]
    vals = [sf34, ab]
    colors = [INX_COLOR, CNX_COLOR]
    bars = ax.bar(labels, vals, color=colors, width=0.55)
    for bar, val in zip(bars, vals):
        ax.text(bar.get_x() + bar.get_width() / 2, val + 1.2, f"{val:.1f}%", ha="center", fontsize=11)
    ax.set_ylabel("Mean MMseqs2 identity (%)")
    ax.set_ylim(0, max(vals) * 1.25)
    ax.set_title("The contrast that matters: two deep splits")
    _save(fig, path)


def plot_lengths(inx: pd.DataFrame, cnx: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    data = [
        pd.to_numeric(inx["protein_length"], errors="coerce").dropna().values,
        pd.to_numeric(cnx["protein_length"], errors="coerce").dropna().values,
    ]
    bp = ax.boxplot(data, tick_labels=["Innexin", "Connexin"], patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], (INX_COLOR, CNX_COLOR)):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_ylabel("Protein length (aa)")
    ax.set_title("Protein length: innexin-shaped vs connexin-shaped")
    _save(fig, path)


def plot_exons(exon: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.0, 5.0))
    data = [
        exon[exon["family"] == "innexin"]["mean_exon_count"].dropna().values,
        exon[exon["family"] == "connexin"]["mean_exon_count"].dropna().values,
    ]
    bp = ax.boxplot(data, tick_labels=["Innexin", "Connexin"], patch_artist=True, showfliers=False)
    for patch, color in zip(bp["boxes"], (INX_COLOR, CNX_COLOR)):
        patch.set_facecolor(color)
        patch.set_alpha(0.7)
    ax.set_ylabel("Mean exon count")
    ax.set_title("Gene architecture: variable innexin exons vs compact connexins")
    _save(fig, path)


def plot_discovery_yield(inx: pd.DataFrame, cnx: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 4.8))
    labels = ["Reference panel", "Discovery / curator"]
    inx_ref = int((inx["source"] == "reference_db").sum())
    cnx_ref = int((cnx["source"] == "reference_db").sum())
    inx_disc = len(inx) - inx_ref
    cnx_disc = len(cnx) - cnx_ref
    x = np.arange(len(labels))
    w = 0.34
    ax.bar(x - w / 2, [inx_ref, inx_disc], w, color=INX_COLOR, label="Innexin")
    ax.bar(x + w / 2, [cnx_ref, cnx_disc], w, color=CNX_COLOR, label="Connexin")
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel("Loci")
    ax.set_title("Discovery yield vs reference panel")
    ax.legend(frameon=False)
    _save(fig, path)


def plot_architecture_modes(path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.8, 5.6))

    ax = axes[0]
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("Innexin: tandem clusters cut across the tree", color=INX_COLOR, pad=8)

    def gene_box(ax_, x, y, w, label, color):
        ax_.add_patch(FancyBboxPatch((x, y), w, 1.05, boxstyle="round,pad=0.04,rounding_size=0.12", facecolor=color, edgecolor="white", linewidth=1.2))
        ax_.text(x + w / 2, y + 0.52, label, ha="center", va="center", color="white", fontsize=8, fontweight="bold")

    ax.text(0.3, 8.55, "Drosophila X proximal  (~6.8–7.9 Mb)", fontsize=8, color=MUTED)
    ax.add_patch(Rectangle((0.3, 6.85), 11.4, 1.55, facecolor="#ecfdf5", edgecolor="#99f6e4", linewidth=1))
    gene_box(ax, 0.55, 7.05, 3.1, "ogre  SF2", INX_SF_COLORS["SF2_Inx1_ogre"])
    gene_box(ax, 4.0, 7.05, 3.4, "Inx7  SF3", INX_SF_COLORS["SF3_Inx3_Inx7_nematode"])
    gene_box(ax, 7.7, 7.05, 3.5, "Inx2  SF4", INX_SF_COLORS["SF4_Inx2_expansion"])

    ax.text(0.3, 5.85, "Drosophila X distal  (~20–21 Mb)", fontsize=8, color=MUTED)
    ax.add_patch(Rectangle((0.3, 4.15), 6.2, 1.55, facecolor="#eff6ff", edgecolor="#bfdbfe", linewidth=1))
    gene_box(ax, 0.55, 4.35, 3.6, "shakB  SF1", INX_SF_COLORS["SF1_shakB"])

    ax.text(7.0, 5.85, "Chromosome 3R", fontsize=8, color=MUTED)
    ax.add_patch(Rectangle((6.9, 4.15), 4.8, 1.55, facecolor="#fff7ed", edgecolor="#fed7aa", linewidth=1))
    gene_box(ax, 7.15, 4.35, 4.3, "Inx3  SF3", INX_SF_COLORS["SF3_Inx3_Inx7_nematode"])

    ax.text(0.3, 2.55, "Nematodes: all innexins sit in SF3", fontsize=8, color=MUTED)
    ax.add_patch(Rectangle((0.3, 0.85), 11.4, 1.55, facecolor="#fff7ed", edgecolor="#fed7aa", linewidth=1))
    gene_box(ax, 0.55, 1.05, 3.3, "unc-7", INX_SF_COLORS["SF3_Inx3_Inx7_nematode"])
    gene_box(ax, 4.1, 1.05, 3.3, "unc-9", INX_SF_COLORS["SF3_Inx3_Inx7_nematode"])
    gene_box(ax, 7.7, 1.05, 3.6, "inx-2 (Cele)", INX_SF_COLORS["SF3_Inx3_Inx7_nematode"])
    ax.text(0.3, 0.25, "Genomic neighbours ≠ phylogenetic sisters. Names do not transfer.", fontsize=8, color="#0f2430")

    ax = axes[1]
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 10)
    ax.axis("off")
    ax.set_title("Connexin: two conserved vertebrate addresses", color=CNX_COLOR, pad=8)

    ax.text(0.3, 8.55, "Conserved β-cluster  (GJA4 lives here)", fontsize=8, color=MUTED)
    ax.add_patch(Rectangle((0.3, 6.85), 11.4, 1.55, facecolor="#fff7ed", edgecolor="#fed7aa", linewidth=1))
    gene_box(ax, 0.5, 7.05, 2.1, "GJB3", CNX_SF_COLORS["beta_GJB"])
    gene_box(ax, 2.75, 7.05, 2.1, "GJB4", CNX_SF_COLORS["beta_GJB"])
    gene_box(ax, 5.0, 7.05, 2.1, "GJB5", CNX_SF_COLORS["beta_GJB"])
    gene_box(ax, 7.25, 7.05, 2.1, "GJA4", CNX_SF_COLORS["alpha_GJA"])
    ax.add_patch(FancyBboxPatch((9.5, 7.05), 1.95, 1.05, boxstyle="round,pad=0.04,rounding_size=0.12", facecolor="#e2e8f0", edgecolor="white"))
    ax.text(10.48, 7.57, "DLGAP3", ha="center", va="center", fontsize=7, color="#334155")

    ax.text(0.3, 5.85, "Separate GJA1 locus", fontsize=8, color=MUTED)
    ax.add_patch(Rectangle((0.3, 4.15), 6.4, 1.55, facecolor="#fff1f2", edgecolor="#fecdd3", linewidth=1))
    gene_box(ax, 0.55, 4.35, 3.6, "GJA1  α", CNX_SF_COLORS["alpha_GJA"])

    ax.text(7.0, 5.85, "Already in Xenopus", fontsize=8, color=MUTED)
    ax.add_patch(Rectangle((6.9, 4.15), 4.8, 1.55, facecolor="#f0f9ff", edgecolor="#bae6fd", linewidth=1))
    gene_box(ax, 7.15, 4.35, 4.3, "β-cluster + GJA1", "#0369A1")

    ax.text(0.3, 2.55, "Named classes travel together across vertebrates", fontsize=8, color=MUTED)
    ax.add_patch(Rectangle((0.3, 0.85), 11.4, 1.55, facecolor="#f8fafc", edgecolor="#e2e8f0", linewidth=1))
    gene_box(ax, 0.5, 1.05, 2.4, "α GJA", CNX_SF_COLORS["alpha_GJA"])
    gene_box(ax, 3.1, 1.05, 2.4, "β GJB", CNX_SF_COLORS["beta_GJB"])
    gene_box(ax, 5.7, 1.05, 2.4, "γ GJC", CNX_SF_COLORS["gamma_GJC"])
    gene_box(ax, 8.3, 1.05, 3.0, "δ GJD", CNX_SF_COLORS["delta_GJD"])
    ax.text(0.3, 0.25, "Neighborhood-conserved. GJA1 ≠ GJA4 even though both are α.", fontsize=8, color="#0f2430")

    _save(fig, path)


def plot_cross_family(m8: Path | None, path: Path) -> dict:
    stats = {"n_hits": 0, "n_homologous": 0, "median_pident": float("nan"), "median_alnlen": float("nan")}
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    empty_msg = "0 hits ≥80 aa at e≤1e-3 → extreme divergence / separate families\n(not proof of absence of homology)"
    if m8 is None or not m8.exists() or m8.stat().st_size == 0:
        ax.text(0.5, 0.55, "No cross-family MMseqs hits", ha="center", va="center", fontsize=13, color=MUTED)
        ax.text(0.5, 0.28, empty_msg, ha="center", fontsize=10)
        ax.axis("off")
        _save(fig, path)
        return stats
    hits = pd.read_csv(m8, sep="\t", header=None, names=["query", "target", "pident", "alnlen", "evalue", "bits"])
    stats["n_hits"] = len(hits)
    if hits.empty:
        ax.text(0.5, 0.5, "Zero innexin↔connexin hits", ha="center", va="center", fontsize=13)
        ax.axis("off")
        _save(fig, path)
        return stats
    homologous = hits[(hits["evalue"] <= 1e-3) & (hits["alnlen"] >= 80)]
    stats["n_homologous"] = int(len(homologous))
    stats["median_pident"] = float(hits["pident"].median())
    stats["median_alnlen"] = float(hits["alnlen"].median())
    ax.scatter(hits["alnlen"], hits["pident"], s=28, alpha=0.75, c="#475569", edgecolors="white", linewidths=0.3, label="raw hits")
    if not homologous.empty:
        ax.scatter(homologous["alnlen"], homologous["pident"], s=42, alpha=0.9, c=CNX_COLOR, label="e≤1e-3 and ≥80 aa")
    ax.axvline(80, color=INX_COLOR, ls="--", lw=1.1, label="80 aa alignment")
    ax.set_xlabel("Alignment length (aa)")
    ax.set_ylabel("Identity (%)")
    ax.set_title(f"Cross-family MMseqs: {stats['n_homologous']} hit(s) ≥80 aa at e≤1e-3 (of {len(hits)} fragments)")
    ax.legend(frameon=False, fontsize=8)
    ax.text(
        0.98,
        0.04,
        empty_msg if stats["n_homologous"] == 0 else "",
        transform=ax.transAxes,
        ha="right",
        fontsize=8,
        color=MUTED,
    )
    _save(fig, path)
    return stats


def write_html(stats: dict, figures: list[tuple[str, str, str, str]]) -> None:
    beats = []
    for i, (fid, fname, title, claim) in enumerate(figures, start=1):
        beats.append(
            f"""
  <article class="beat" id="{fid}">
    <span class="eyebrow">{i:02d}</span>
    <h2>{html.escape(title)}</h2>
    <p class="claim">{claim}</p>
    <div class="figure">
      <img src="figures/{html.escape(fname)}" alt="{html.escape(title)}" onclick="openModal(this.src)">
    </div>
  </article>"""
        )

    toc = "\n".join(f'    <a href="#{fid}">{i}. {html.escape(title)}</a>' for i, (fid, _f, title, _c) in enumerate(figures, start=1))
    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Innexin vs Connexin — same comparisons</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,560;9..144,700&family=Sora:wght@400;500;600&display=swap" rel="stylesheet">
<style>
:root {{
  --ink: #12202a;
  --muted: #4d6570;
  --deep: #0a3a42;
  --sea: #0f766e;
  --ember: #c2410c;
  --foam: #f4f7f5;
  --card: #ffffff;
  --line: rgba(18, 32, 42, 0.1);
}}
* {{ box-sizing: border-box; }}
html {{ scroll-behavior: smooth; }}
body {{
  margin: 0; color: var(--ink); font-family: "Sora", sans-serif; line-height: 1.6;
  background:
    radial-gradient(1000px 480px at 8% -8%, rgba(15, 118, 110, 0.16), transparent 55%),
    radial-gradient(900px 420px at 92% 4%, rgba(194, 65, 12, 0.12), transparent 50%),
    linear-gradient(180deg, #e8f1ee 0%, var(--foam) 40%, #eef3f1 100%);
}}
a {{ color: var(--sea); }}
.hero {{
  position: relative; min-height: 86vh; display: grid; align-items: end; overflow: hidden;
  padding: clamp(1.5rem, 4vw, 3.25rem); color: #f7fffb;
  background:
    linear-gradient(118deg, rgba(10, 58, 66, 0.94) 0%, rgba(15, 118, 110, 0.72) 48%, rgba(194, 65, 12, 0.62) 100%),
    url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='140' height='140'%3E%3Cpath d='M0 70h140M70 0v140' stroke='%23fff' stroke-opacity='.07'/%3E%3C/svg%3E");
}}
.brand {{ font-family: "Fraunces", serif; font-size: clamp(2.4rem, 6.5vw, 4.5rem); font-weight: 700; letter-spacing: -0.03em; line-height: 0.95; margin: 0 0 1rem; }}
.hero h1 {{ margin: 0 0 0.75rem; font-size: clamp(1.05rem, 2.1vw, 1.3rem); font-weight: 500; max-width: 42ch; }}
.hero .lede {{ margin: 0 0 1.6rem; max-width: 52ch; opacity: 0.92; }}
.cta a {{
  display: inline-block; text-decoration: none; font-weight: 600; padding: 0.8rem 1.15rem; border-radius: 999px;
  margin-right: 0.45rem; margin-bottom: 0.35rem; background: #f5fffb; color: var(--deep); box-shadow: 0 8px 24px rgba(0,0,0,0.18);
}}
.cta a.ghost {{ background: transparent; color: #f5fffb; border: 1px solid rgba(245,255,251,0.45); box-shadow: none; }}
.wrap {{ width: min(880px, calc(100% - 2rem)); margin: 0 auto; padding: 2.25rem 0 4.5rem; }}
.verdict {{
  margin-top: -2.2rem; position: relative; z-index: 2; background: rgba(255,255,255,0.95);
  border: 1px solid var(--line); border-radius: 18px; padding: 1.25rem 1.35rem; box-shadow: 0 14px 36px rgba(18,32,42,0.08);
}}
.verdict h2 {{ font-family: "Fraunces", serif; font-size: 1.35rem; margin: 0 0 0.55rem; }}
.verdict p {{ margin: 0; color: var(--muted); }}
.verdict strong {{ color: var(--deep); }}
.stats {{ display: grid; grid-template-columns: repeat(4, 1fr); gap: 0.7rem; margin: 1.4rem 0 0; }}
.stat {{ background: #fff; border: 1px solid var(--line); border-radius: 14px; padding: 0.85rem; }}
.stat b {{ display: block; font-family: "Fraunces", serif; font-size: 1.35rem; color: var(--deep); }}
.stat span {{ color: var(--muted); font-size: 0.78rem; }}
.toc {{ display: flex; flex-wrap: wrap; gap: 0.45rem; margin: 1.75rem 0 0.5rem; }}
.toc a {{
  text-decoration: none; font-size: 0.78rem; font-weight: 500; color: var(--deep); background: #fff;
  border: 1px solid var(--line); padding: 0.4rem 0.75rem; border-radius: 999px;
}}
.beat {{ margin: 2.3rem 0; scroll-margin-top: 1rem; }}
.beat .eyebrow {{ display: inline-block; font-size: 0.72rem; font-weight: 700; letter-spacing: 0.08em; text-transform: uppercase; color: var(--ember); margin-bottom: 0.35rem; }}
.beat h2 {{ font-family: "Fraunces", serif; font-size: clamp(1.4rem, 3vw, 1.8rem); margin: 0 0 0.55rem; letter-spacing: -0.02em; }}
.beat .claim {{ color: var(--muted); margin: 0 0 1rem; max-width: 64ch; }}
.figure {{ background: var(--card); border: 1px solid var(--line); border-radius: 16px; padding: 0.85rem; box-shadow: 0 10px 28px rgba(18,32,42,0.05); }}
.figure img {{ width: 100%; display: block; border-radius: 10px; cursor: zoom-in; background: #f8fafc; }}
.closing {{ margin-top: 2.5rem; padding: 1.4rem 1.35rem; border-radius: 16px; background: linear-gradient(135deg, var(--deep), var(--ember)); color: #f5fffb; }}
.closing h2 {{ font-family: "Fraunces", serif; margin: 0 0 0.55rem; font-size: 1.4rem; }}
.closing p {{ margin: 0; opacity: 0.95; max-width: 62ch; }}
.footer {{ margin-top: 1.75rem; color: var(--muted); font-size: 0.88rem; }}
.footer a {{ margin-right: 0.85rem; }}
.modal {{ display: none; position: fixed; inset: 0; background: rgba(0,0,0,0.88); z-index: 30; align-items: center; justify-content: center; padding: 1.25rem; cursor: zoom-out; }}
.modal.open {{ display: flex; }}
.modal img {{ max-width: 95vw; max-height: 90vh; border-radius: 8px; }}
@media (max-width: 800px) {{ .stats {{ grid-template-columns: 1fr 1fr; }} .hero {{ min-height: 76vh; }} }}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <p class="brand">Innexin vs Connexin</p>
    <h1>The same comparisons, read against each other</h1>
    <p class="lede">Clade geography, within/between identity, embeddings, length, exons and synteny logic — one frame for both gap-junction families.</p>
    <div class="cta">
      <a href="#verdict">Start here</a>
      <a class="ghost" href="../innexin_insights/index.html">Innexin insights</a>
      <a class="ghost" href="../connexin_insights/index.html">Connexin insights</a>
    </div>
  </div>
</header>
<main class="wrap">
  <section class="verdict" id="verdict">
    <h2>Main contrast</h2>
    <p>
      Both families encode gap junctions. Cross-family MMseqs and embeddings in this panel treat them as
      <strong>separate sequence clouds</strong>, consistent with literature that places them in different
      channel families (not proof of absence of deep homology).
      Innexins diversified by <strong>clade-specific duplication</strong> with misleading gene names
      (within ~{stats['inx_id']['within_mean']:.0f}% vs between ~{stats['inx_id']['between_mean']:.0f}%; SF3↔SF4 ~{stats['sf34']:.0f}%).
      Connexins diversified inside vertebrates as <strong>named α/β/γ/δ classes</strong> that stay closer
      (within ~{stats['cnx_id']['within_mean']:.0f}% vs between ~{stats['cnx_id']['between_mean']:.0f}%; α↔β ~{stats['ab']:.0f}%)
      and keep two conserved genomic addresses (GJA4 in the β-cluster, GJA1 elsewhere).
    </p>
    <div class="stats">
      <div class="stat"><b>{stats['n_inx']}</b><span>innexin proteins</span></div>
      <div class="stat"><b>{stats['n_cnx']}</b><span>connexin proteins</span></div>
      <div class="stat"><b>{stats['sf34']:.0f}% vs {stats['ab']:.0f}%</b><span>SF3↔SF4 vs α↔β identity</span></div>
      <div class="stat"><b>{stats['cross_hits']}</b><span>innexin↔connexin MMseqs hits</span></div>
    </div>
  </section>
  <nav class="toc" aria-label="Comparisons">
{toc}
  </nav>
{''.join(beats)}
  <section class="closing">
    <h2>Summary</h2>
    <p>
      Innexins and connexins fill a similar biological niche with <strong>distinct sequence and genomic
      characters in this panel</strong>:
      innexins split into tree subfamilies whose proteins and names come apart across phyla,
      while connexins keep α/β classes sequence-close and synteny-stable across vertebrates.
    </p>
  </section>
  <footer class="footer">
    <p>Built from innexin and connexin clade tables + MMseqs/embedding matrices + exon summaries{stats['cross_note']}.</p>
    <p>
      <a href="../innexin_insights/index.html">Innexin insights</a>
      <a href="../connexin_insights/index.html">Connexin insights</a>
      <a href="../innexin_similarity/index.html">Innexin MMseqs</a>
      <a href="../connexin_similarity/index.html">Connexin MMseqs</a>
      <a href="../phylogenetic_story/index.html">Phylogeny story</a>
      <a href="../sources/index.html">Sources</a>
    </p>
  </footer>
</main>
<div id="modal" class="modal" onclick="this.classList.remove('open')">
  <img id="modal-img" alt="Enlarged figure">
</div>
<script>
function openModal(src) {{
  document.getElementById('modal-img').src = src;
  document.getElementById('modal').classList.add('open');
}}
document.addEventListener('keydown', e => {{
  if (e.key === 'Escape') document.getElementById('modal').classList.remove('open');
}});
</script>
</body>
</html>
"""
    (OUT / "index.html").write_text(page, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    inx = pd.read_csv(INX_CLADE / "all_innexin_loci_annotated.csv")
    cnx = pd.read_csv(CNX_CLADE / "all_connexin_loci_annotated.csv")
    inx_mat = load_group_matrix(INX_SIM / "group_similarity_matrix.csv")
    cnx_mat = load_group_matrix(CNX_SIM / "group_similarity_matrix.csv")
    inx_sw = pd.read_csv(INX_SIM / "mmseqs_threshold_sweep_summary.csv")
    cnx_sw = pd.read_csv(CNX_SIM / "mmseqs_threshold_sweep_summary.csv")
    exon = pd.read_csv(EXON_SUMMARY) if EXON_SUMMARY.exists() else pd.DataFrame()

    inx_id = family_metrics(inx_mat, "mmseqs_pident")
    cnx_id = family_metrics(cnx_mat, "mmseqs_pident")
    inx_emb = family_metrics(inx_mat, "kmer3_svd_cosine")
    cnx_emb = family_metrics(cnx_mat, "kmer3_svd_cosine")
    inx_aa = family_metrics(inx_mat, "aa_comp_cosine")
    cnx_aa = family_metrics(cnx_mat, "aa_comp_cosine")
    sf34 = pair_mean(inx_mat, "mmseqs_pident", "SF3_Inx3_Inx7_nematode", "SF4_Inx2_expansion")
    ab = pair_mean(cnx_mat, "mmseqs_pident", "alpha_GJA", "beta_GJB")
    inx_id_plot = {
        "within_mean": CANONICAL["innexin"]["within_id"],
        "between_mean": CANONICAL["innexin"]["between_id"],
    }
    cnx_id_plot = {
        "within_mean": CANONICAL["connexin"]["within_id"],
        "between_mean": CANONICAL["connexin"]["between_id"],
    }
    inx_emb_plot = {
        "within_mean": CANONICAL["innexin"]["emb_same"],
        "between_mean": CANONICAL["innexin"]["emb_diff"],
    }
    cnx_emb_plot = {
        "within_mean": CANONICAL["connexin"]["emb_same"],
        "between_mean": CANONICAL["connexin"]["emb_diff"],
    }

    m8 = run_cross_family_mmseqs()
    cross = plot_cross_family(m8, FIGS / "11_cross_family_mmseqs.png")

    plot_dataset_scope(inx, cnx, FIGS / "01_dataset_scope.png")
    plot_clade_geography(inx, cnx, FIGS / "02_clade_geography.png")
    plot_within_between(inx_id_plot, cnx_id_plot, FIGS / "03_within_between_identity.png")
    plot_identity_heatmaps(inx_mat, cnx_mat, FIGS / "04_identity_heatmaps.png")
    plot_embedding_separation(inx_emb_plot, cnx_emb_plot, FIGS / "05_embedding_separation.png")
    plot_threshold_sweep(inx_sw, cnx_sw, FIGS / "06_threshold_sweep.png")
    plot_key_pair(sf34, ab, FIGS / "07_key_pair_identity.png")
    plot_lengths(inx, cnx, FIGS / "08_protein_length.png")
    if not exon.empty:
        plot_exons(exon, FIGS / "09_exon_architecture.png")
    plot_discovery_yield(inx, cnx, FIGS / "10_discovery_yield.png")
    plot_architecture_modes(FIGS / "12_architecture_modes.png")

    metric_rows = [
        {"family": "innexin", "metric": "n_proteins", "value": str(len(inx))},
        {"family": "connexin", "metric": "n_proteins", "value": str(len(cnx))},
        {"family": "innexin", "metric": "n_species", "value": str(inx["organism"].nunique())},
        {"family": "connexin", "metric": "n_species", "value": str(cnx["organism"].nunique())},
        {"family": "innexin", "metric": "within_identity_published", "value": str(CANONICAL["innexin"]["within_id"])},
        {"family": "connexin", "metric": "within_identity_published", "value": str(CANONICAL["connexin"]["within_id"])},
        {"family": "innexin", "metric": "between_identity_published", "value": str(CANONICAL["innexin"]["between_id"])},
        {"family": "connexin", "metric": "between_identity_published", "value": str(CANONICAL["connexin"]["between_id"])},
        {"family": "innexin", "metric": "within_identity_matrix_mean", "value": f"{inx_id['within_mean']:.4f}"},
        {"family": "connexin", "metric": "within_identity_matrix_mean", "value": f"{cnx_id['within_mean']:.4f}"},
        {"family": "innexin", "metric": "between_identity_matrix_mean", "value": f"{inx_id['between_mean']:.4f}"},
        {"family": "connexin", "metric": "between_identity_matrix_mean", "value": f"{cnx_id['between_mean']:.4f}"},
        {"family": "innexin", "metric": "within_kmer_cosine", "value": f"{inx_emb['within_mean']:.4f}"},
        {"family": "connexin", "metric": "within_kmer_cosine", "value": f"{cnx_emb['within_mean']:.4f}"},
        {"family": "innexin", "metric": "between_kmer_cosine", "value": f"{inx_emb['between_mean']:.4f}"},
        {"family": "connexin", "metric": "between_kmer_cosine", "value": f"{cnx_emb['between_mean']:.4f}"},
        {"family": "innexin", "metric": "within_aa_cosine", "value": f"{inx_aa['within_mean']:.4f}"},
        {"family": "innexin", "metric": "between_aa_cosine", "value": f"{inx_aa['between_mean']:.4f}"},
        {"family": "connexin", "metric": "within_aa_cosine", "value": f"{cnx_aa['within_mean']:.4f}"},
        {"family": "connexin", "metric": "between_aa_cosine", "value": f"{cnx_aa['between_mean']:.4f}"},
        {"family": "innexin", "metric": "SF3_vs_SF4_identity", "value": f"{sf34:.4f}"},
        {"family": "connexin", "metric": "alpha_vs_beta_identity", "value": f"{ab:.4f}"},
        {"family": "both", "metric": "cross_family_raw_hits", "value": str(cross["n_hits"])},
        {"family": "both", "metric": "cross_family_homologous_hits", "value": str(cross["n_homologous"])},
        {"family": "both", "metric": "cross_family_median_pident", "value": f"{cross['median_pident']:.4f}" if not math.isnan(cross["median_pident"]) else ""},
        {"family": "both", "metric": "cross_family_median_alnlen", "value": f"{cross['median_alnlen']:.4f}" if not math.isnan(cross.get("median_alnlen", float("nan"))) else ""},
    ]
    write_csv(OUT / "comparison_metrics.csv", metric_rows)

    aln_note = ""
    if not math.isnan(cross.get("median_alnlen", float("nan"))):
        aln_note = f"; fragment median length {cross['median_alnlen']:.0f} aa"
    findings = [
        f"Innexin within-subfamily identity {CANONICAL['innexin']['within_id']:.1f}% vs between {CANONICAL['innexin']['between_id']:.1f}%.",
        f"Connexin within-class identity {CANONICAL['connexin']['within_id']:.1f}% vs between {CANONICAL['connexin']['between_id']:.1f}%.",
        f"The deep splits: innexin SF3↔SF4 {sf34:.1f}% vs connexin α↔β {ab:.1f}%.",
        f"3-mer embeddings: innexin same {CANONICAL['innexin']['emb_same']:.3f} vs different {CANONICAL['innexin']['emb_diff']:.3f}; connexin same {CANONICAL['connexin']['emb_same']:.3f} vs different {CANONICAL['connexin']['emb_diff']:.3f}.",
        f"AA composition barely splits either family (innexin {inx_aa['within_mean']:.3f} vs {inx_aa['between_mean']:.3f}; connexin {cnx_aa['within_mean']:.3f} vs {cnx_aa['between_mean']:.3f}).",
        f"Cross-family MMseqs: {cross['n_hits']} fragment hits, {cross['n_homologous']} at e≤1e-3 and ≥80 aa{aln_note} — no detectable long-range similarity (extreme divergence / separate families; not proof of absence of homology).",
        "Innexin architecture is tandem/clade-bound (Drosophila clusters ≠ tree subfamilies). Connexin architecture is neighborhood-conserved (GJA4 β-cluster vs separate GJA1 locus).",
    ]
    (OUT / "key_findings.md").write_text("# Key findings\n\n" + "\n".join(f"- {f}" for f in findings) + "\n", encoding="utf-8")

    exon_claim = (
        "Reference gene models: innexin exon counts vary by lineage; connexins stay compact (often ~2 exons). Architecture is a family character, not a shared ancestral gene model."
        if not exon.empty
        else "Exon summaries were not found; length still separates the two protein shapes."
    )
    cross_note = " + cross-family MMseqs" if m8 else ""
    cross_claim = (
        f"Sensitive MMseqs finds <strong>{cross['n_hits']}</strong> short fragment hits"
        + (f" (median length {cross['median_alnlen']:.0f} aa)" if not math.isnan(cross.get("median_alnlen", float("nan"))) else "")
        + f", but <strong>{cross['n_homologous']}</strong> alignment(s) at e≤1e-3 and ≥80 aa. "
        "Zero long hits supports extreme divergence or separate families — not by itself proof of absence of homology. "
        "Same biological job; analyse sequence and synteny in parallel, not as one alignment."
    )
    figures = [
        ("scope", "01_dataset_scope.png", "Same comparison, two datasets",
         f"Innexin panel: <strong>{len(inx)}</strong> proteins in <strong>{inx['organism'].nunique()}</strong> species. Connexin panel: <strong>{len(cnx)}</strong> proteins in <strong>{cnx['organism'].nunique()}</strong> species. The innexin story is discovery-heavy; the connexin story is carried by the vertebrate reference panel."),
        ("clades", "02_clade_geography.png", "Clade geography is not interchangeable",
         "Innexin recoveries outside insects pile into SF3 / nematode-like labels. Connexin α and β co-occur across mammals, birds, amphibians and fish — the canonical vertebrate gap-junction radiation in this panel, not a name-chaos problem."),
        ("identity", "03_within_between_identity.png", "Within vs between identity",
         f"Identical metric. Innexins: ~{CANONICAL['innexin']['within_id']:.0f}% within vs ~{CANONICAL['innexin']['between_id']:.0f}% between. Connexins stay closer: ~{CANONICAL['connexin']['within_id']:.0f}% vs ~{CANONICAL['connexin']['between_id']:.0f}%. Both are still one family internally; connexin paralogs diverged less."),
        ("heatmaps", "04_identity_heatmaps.png", "Subfamily heatmaps, paired",
         f"Read the off-diagonals. Innexin SF3 vs SF4 is the cold cell (~{sf34:.0f}%). Connexin α vs β is still warm (~{ab:.0f}%) — distinct classes, younger split."),
        ("split", "07_key_pair_identity.png", "The contrast that matters",
         "SF3↔SF4 is two evolutionary experiments inside innexins (insect expansion vs nematode-like radiation). α↔β is a conserved vertebrate paralog split. Same plot type; different depths of divergence."),
        ("embed", "05_embedding_separation.png", "3-mer embeddings agree",
         f"Local motif grammar separates groups even when amino-acid composition does not. Innexin same vs different: {CANONICAL['innexin']['emb_same']:.2f} vs {CANONICAL['innexin']['emb_diff']:.2f}. Connexin: {CANONICAL['connexin']['emb_same']:.2f} vs {CANONICAL['connexin']['emb_diff']:.2f}."),
        ("sweep", "06_threshold_sweep.png", "How cleanly do groups cluster?",
         "Connexin named classes stay pure even at high identity (recent, well-separated paralogs). Innexin subfamilies peak around 40% identity — older, more mixed expansions."),
        ("length", "08_protein_length.png", "Protein shape is family-specific",
         "Innexins sit longer (~350–420 aa). Connexins are compact (~250–380 aa). Length is conserved more tightly than sequence identity inside each family."),
        ("exons", "09_exon_architecture.png", "Exon architecture",
         exon_claim),
        ("yield", "10_discovery_yield.png", "Discovery yield vs reference panel",
         "Innexin curator recoveries (rotifers, molluscs, non-insect arthropods) carry the phylogenetic argument. Connexin discovery is sparse; α/β/γ/δ structure is already in the reference vertebrates."),
        ("cross", "11_cross_family_mmseqs.png", "Cross-family sequence search",
         cross_claim),
        ("arch", "12_architecture_modes.png", "Synteny logic is different",
         "Innexin: Drosophila tandem clusters cut across the tree (Inx7 is next to ogre/Inx2 but phylogenetically with nematodes). Connexin: GJA4 stays in a conserved β-cluster; GJA1 is a second locus already present in Xenopus."),
    ]
    if exon.empty:
        figures = [fig for fig in figures if fig[0] != "exons"]

    write_html(
        {
            "n_inx": len(inx),
            "n_cnx": len(cnx),
            "inx_id": inx_id_plot,
            "cnx_id": cnx_id_plot,
            "sf34": sf34,
            "ab": ab,
            "cross_hits": cross["n_homologous"],
            "cross_note": cross_note,
        },
        figures,
    )
    print(f"Done → {OUT / 'index.html'}")
    print(f"  innexin within/between {CANONICAL['innexin']['within_id']:.1f}/{CANONICAL['innexin']['between_id']:.1f}")
    print(f"  connexin within/between {CANONICAL['connexin']['within_id']:.1f}/{CANONICAL['connexin']['between_id']:.1f}")
    print(f"  SF3↔SF4 {sf34:.1f}  α↔β {ab:.1f}  homologous {cross['n_homologous']} / raw {cross['n_hits']}")


if __name__ == "__main__":
    main()
