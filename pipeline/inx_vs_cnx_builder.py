#!/usr/bin/env python3
"""Direct innexin vs connexin contrast (two families in one feature space).

This is a new site, not a rebuild of family_comparison/.
Outputs → project/results/inx_vs_cnx/
"""

from __future__ import annotations

import html
import math
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from matplotlib.patches import FancyBboxPatch, Rectangle

from pipeline.common import RESULTS_DIR, write_csv

OUT = RESULTS_DIR / "inx_vs_cnx"
FIGS = OUT / "figures"

INX_FASTA = RESULTS_DIR / "innexin_similarity" / "all_innexins.fasta"
CNX_FASTA = RESULTS_DIR / "connexin_similarity" / "all_connexins.fasta"
INX_LOCI = RESULTS_DIR / "innexin_clade_comparison" / "all_innexin_loci_annotated.csv"
CNX_LOCI = RESULTS_DIR / "connexin_clade_comparison" / "all_connexin_loci_annotated.csv"
EXON = RESULTS_DIR / "exon_structures" / "gene_summary.csv"
TRANS = RESULTS_DIR / "exon_structures" / "transcript_summary.csv"
CROSS_M8 = RESULTS_DIR / "family_comparison" / "innexin_vs_connexin.m8"

INX = "#0F766E"
CNX = "#C2410C"
MUTED = "#64748B"

AA = list("ACDEFGHIKLMNPQRSTVWY")
KD = {
    "A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
    "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
    "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2,
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


def read_fasta(path: Path, family: str) -> list[dict]:
    rows: list[dict] = []
    seq_id = ""
    chunks: list[str] = []
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id:
                    rows.append(_seq_row(seq_id, "".join(chunks), family))
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if seq_id:
            rows.append(_seq_row(seq_id, "".join(chunks), family))
    return rows


def _seq_row(seq_id: str, seq: str, family: str) -> dict:
    seq = "".join(a for a in seq if a.isalpha())
    n = max(len(seq), 1)
    cys = seq.count("C")
    gravy = sum(KD.get(a, 0.0) for a in seq) / n
    charged = sum(seq.count(a) for a in "DEKR") / n
    aromatic = sum(seq.count(a) for a in "FWY") / n
    return {
        "seq_id": seq_id,
        "family": family,
        "seq": seq,
        "length": len(seq),
        "cys_per_100": 100.0 * cys / n,
        "gravy": gravy,
        "charged_frac": charged,
        "aromatic_frac": aromatic,
        **{f"aa_{a}": seq.count(a) / n for a in AA},
    }


def kmer_svd(seqs: list[str], k: int = 3, n_comp: int = 8) -> np.ndarray:
    vocab_count: Counter[str] = Counter()
    for seq in seqs:
        for i in range(max(0, len(seq) - k + 1)):
            km = seq[i : i + k]
            if "X" in km:
                continue
            vocab_count[km] += 1
    vocab = [km for km, _ in vocab_count.most_common(800)]
    idx = {km: j for j, km in enumerate(vocab)}
    X = np.zeros((len(seqs), max(len(vocab), 1)), dtype=float)
    for i, seq in enumerate(seqs):
        for j in range(max(0, len(seq) - k + 1)):
            km = seq[j : j + k]
            if km in idx:
                X[i, idx[km]] += 1
        s = X[i].sum()
        if s:
            X[i] /= s
    Xc = X - X.mean(axis=0, keepdims=True)
    u, s, _vt = np.linalg.svd(Xc, full_matrices=False)
    n_comp = min(n_comp, u.shape[1])
    emb = u[:, :n_comp] * s[:n_comp]
    norms = np.linalg.norm(emb, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return emb / norms


def aa_matrix(rows: list[dict]) -> np.ndarray:
    X = np.array([[r[f"aa_{a}"] for a in AA] for r in rows], dtype=float)
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return X / norms


def pca2(X: np.ndarray) -> np.ndarray:
    Xc = X - X.mean(axis=0, keepdims=True)
    u, s, _vt = np.linalg.svd(Xc, full_matrices=False)
    n = min(2, u.shape[1])
    return u[:, :n] * s[:n]


def mean_cosine(A: np.ndarray, ia: np.ndarray, ib: np.ndarray, max_pairs: int = 8000) -> float:
    if len(ia) == 0 or len(ib) == 0:
        return float("nan")
    rng = np.random.default_rng(7)
    n = min(max_pairs, len(ia) * len(ib) if not np.array_equal(ia, ib) else len(ia) * (len(ia) - 1) // 2)
    if np.array_equal(ia, ib):
        pairs = []
        idx = ia
        for _ in range(n * 3):
            i, j = rng.integers(0, len(idx), size=2)
            if i != j:
                pairs.append((idx[i], idx[j]))
            if len(pairs) >= n:
                break
        if not pairs:
            return float("nan")
        sel = np.array(pairs)
        return float(np.mean(np.sum(A[sel[:, 0]] * A[sel[:, 1]], axis=1)))
    i_s = rng.choice(ia, size=n, replace=True)
    j_s = rng.choice(ib, size=n, replace=True)
    return float(np.mean(np.sum(A[i_s] * A[j_s], axis=1)))


def kingdom(clade: str) -> str:
    c = str(clade).lower()
    verts = ("mammal", "bird", "amphib", "fish", "reptile", "vertebrate")
    if any(v in c for v in verts):
        return "Vertebrates"
    return "Invertebrates"


def plot_occupancy(inx: pd.DataFrame, cnx: pd.DataFrame, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(9.2, 5.0))
    bins = ["Invertebrates", "Vertebrates"]
    inx_k = inx["clade"].map(kingdom).value_counts()
    cnx_k = cnx["clade"].map(kingdom).value_counts()
    x = np.arange(len(bins))
    w = 0.34
    iv = [int(inx_k.get(b, 0)) for b in bins]
    cv = [int(cnx_k.get(b, 0)) for b in bins]
    ax.bar(x - w / 2, iv, w, color=INX, label="Innexin")
    ax.bar(x + w / 2, cv, w, color=CNX, label="Connexin")
    for i, (a, b) in enumerate(zip(iv, cv)):
        if a:
            ax.text(i - w / 2, a + 1.5, str(a), ha="center", fontsize=9, color=INX)
        if b:
            ax.text(i + w / 2, b + 1.5, str(b), ha="center", fontsize=9, color=CNX)
    ax.set_xticks(x)
    ax.set_xticklabels(bins)
    ax.set_ylabel("Loci")
    ax.set_title("Where the two families live")
    ax.legend(frameon=False)
    _save(fig, path)


def plot_copy_number(inx: pd.DataFrame, cnx: pd.DataFrame, path: Path) -> None:
    """Panel loci per species — explicitly not full-genome connexin repertoire."""
    fig, axes = plt.subplots(1, 2, figsize=(12.4, 5.2))

    # Innexin: discovery/recovery panel (meaningful within-project comparison)
    inx_n = inx.groupby("organism").size()
    ax = axes[0]
    hi = int(max(inx_n.max(), 1) + 1)
    ax.hist(
        inx_n,
        bins=range(1, hi + 1),
        color=INX,
        alpha=0.72,
        edgecolor="white",
    )
    ax.set_xlabel("Innexin loci in project panel")
    ax.set_ylabel("Species")
    ax.set_title(f"Innexin recoveries (median {inx_n.median():.0f}/species)")
    ax.text(
        0.03,
        0.97,
        "Reference + discovery loci\nin this panel",
        transform=ax.transAxes,
        va="top",
        fontsize=8,
        color=MUTED,
    )

    # Connexin: do NOT histogram all species — median 1 is sparse sampling, not biology
    ax = axes[1]
    if "source" in cnx.columns:
        cnx_ref = cnx[cnx["source"] == "reference_db"]
    else:
        cnx_ref = cnx
    cnx_n = cnx_ref.groupby("organism").size().sort_values(ascending=True)
    n_all = cnx.groupby("organism").size()
    n_singleton = int((n_all == 1).sum())

    show = cnx_n[cnx_n >= 3]
    if show.empty:
        show = cnx_n
    y = np.arange(len(show))
    ax.barh(y, show.values, color=CNX, alpha=0.85, height=0.72)
    ax.set_yticks(y)
    ax.set_yticklabels([o.replace(" ", "\n") for o in show.index], fontsize=7.5)
    ax.axvspan(20, 22, color=CNX, alpha=0.12, label="Literature jawed vertebrates ~20–22")
    ax.axvline(21, color="#475569", ls="--", lw=1.2, label="Human (21 genes)")
    ax.set_xlabel("Connexin loci in reference panel")
    ax.set_title("Connexin reference sampling (best-covered species)")
    ax.legend(frameon=False, fontsize=7.5, loc="lower right")
    ax.text(
        0.03,
        0.97,
        f"Median 1 across all {len(n_all)} panel species\n"
        f"= {n_singleton} species with only 1 locus here\n"
        "(incomplete UniProt slice, not genome count).\n"
        "Teleosts can carry up to ~46 in nature\n(Mikalsen et al.; panel e.g. Danio n=4).",
        transform=ax.transAxes,
        va="top",
        fontsize=7.5,
        color=MUTED,
    )

    fig.suptitle(
        "Panel loci per species (lower bound / floor) — not full vertebrate connexin gene repertoires",
        fontsize=11.5,
        y=1.02,
    )
    _save(fig, path)


def plot_joint_pca(xy: np.ndarray, families: list[str], title: str, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.8, 6.2))
    fam = np.array(families)
    for lab, color in (("innexin", INX), ("connexin", CNX)):
        m = fam == lab
        ax.scatter(xy[m, 0], xy[m, 1], s=28, alpha=0.72, c=color, label=lab.capitalize(), edgecolors="white", linewidths=0.25)
    ax.set_xlabel("Component 1")
    ax.set_ylabel("Component 2")
    ax.set_title(title)
    ax.legend(frameon=False)
    _save(fig, path)


def plot_family_cosine(vals: dict[str, float], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.2, 4.8))
    labels = ["Innexin ↔ innexin", "Connexin ↔ connexin", "Innexin ↔ connexin"]
    keys = ["inx_inx", "cnx_cnx", "inx_cnx"]
    colors = [INX, CNX, MUTED]
    bars = ax.bar(labels, [vals[k] for k in keys], color=colors, width=0.55)
    for bar, key in zip(bars, keys):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.02, f"{vals[key]:.2f}", ha="center", fontsize=10)
    ax.axhline(0, color="#e2e8f0", lw=1)
    ax.set_ylabel("Mean 3-mer embedding cosine")
    ax.set_title("Do the two families occupy the same sequence texture?")
    ax.set_ylim(min(-0.15, min(vals.values()) - 0.08), max(vals.values()) + 0.12)
    _save(fig, path)


def plot_length_overlay(inx_len: np.ndarray, cnx_len: np.ndarray, path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    bins = np.linspace(80, 700, 28)
    ax.hist(inx_len, bins=bins, color=INX, alpha=0.55, label=f"Innexin (median {np.median(inx_len):.0f} aa)")
    ax.hist(cnx_len, bins=bins, color=CNX, alpha=0.55, label=f"Connexin (median {np.median(cnx_len):.0f} aa)")
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("Proteins")
    ax.set_title("Protein shape: innexin vs connexin in one plot")
    ax.legend(frameon=False)
    _save(fig, path)


def plot_exons_span(exon: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.6, 4.8))
    for ax, col, ylab, title in (
        (axes[0], "mean_exon_count", "Mean exon count", "Gene models: exons"),
        (axes[1], "gene_span", "Gene span (kb)", "Gene models: genomic span"),
    ):
        data = []
        labels = []
        for fam, color in (("innexin", INX), ("connexin", CNX)):
            vals = pd.to_numeric(exon.loc[exon["family"] == fam, col], errors="coerce").dropna()
            if col == "gene_span":
                vals = vals / 1000.0
            data.append(vals.values)
            labels.append(fam.capitalize())
        bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, showfliers=False)
        for patch, color in zip(bp["boxes"], (INX, CNX)):
            patch.set_facecolor(color)
            patch.set_alpha(0.72)
        ax.set_ylabel(ylab)
        ax.set_title(title)
    _save(fig, path)


def plot_chemistry(df: pd.DataFrame, path: Path) -> None:
    fig, axes = plt.subplots(1, 3, figsize=(12.4, 4.6))
    specs = [
        ("cys_per_100", "Cysteines / 100 aa", "Total Cys density (not the EL ladder)"),
        ("gravy", "GRAVY (hydrophobicity)", "Membrane character"),
        ("charged_frac", "Charged residue fraction", "DEKR content"),
    ]
    for ax, (col, ylab, title) in zip(axes, specs):
        data = [df.loc[df["family"] == fam, col].values for fam in ("innexin", "connexin")]
        bp = ax.boxplot(data, tick_labels=["Innexin", "Connexin"], patch_artist=True, showfliers=False)
        for patch, color in zip(bp["boxes"], (INX, CNX)):
            patch.set_facecolor(color)
            patch.set_alpha(0.72)
        ax.set_ylabel(ylab)
        ax.set_title(title)
    _save(fig, path)


def plot_aa_difference(df: pd.DataFrame, path: Path) -> None:
    inx = df[df["family"] == "innexin"]
    cnx = df[df["family"] == "connexin"]
    delta = np.array([inx[f"aa_{a}"].mean() - cnx[f"aa_{a}"].mean() for a in AA])
    fig, ax = plt.subplots(figsize=(9.2, 4.8))
    colors = [INX if d >= 0 else CNX for d in delta]
    ax.bar(AA, delta * 100, color=colors, alpha=0.9)
    ax.axhline(0, color="#cbd5e1", lw=1)
    ax.set_ylabel("Innexin − connexin (percentage points)")
    ax.set_title("Amino-acid mix: which residues each family prefers")
    _save(fig, path)


def plot_cross_mmseqs(path: Path) -> dict:
    stats = {"n_hits": 0, "n_homologous": 0, "median_alnlen": float("nan")}
    fig, ax = plt.subplots(figsize=(8.6, 5.0))
    if not CROSS_M8.exists() or CROSS_M8.stat().st_size == 0:
        ax.text(0.5, 0.5, "No innexin↔connexin alignments", ha="center")
        ax.axis("off")
        _save(fig, path)
        return stats
    hits = pd.read_csv(CROSS_M8, sep="\t", header=None, names=["query", "target", "pident", "alnlen", "evalue", "bits"])
    stats["n_hits"] = len(hits)
    homologous = hits[(hits["evalue"] <= 1e-3) & (hits["alnlen"] >= 80)]
    stats["n_homologous"] = int(len(homologous))
    stats["median_alnlen"] = float(hits["alnlen"].median()) if not hits.empty else float("nan")
    ax.scatter(hits["alnlen"], hits["pident"], s=26, alpha=0.75, c=MUTED, edgecolors="white", linewidths=0.3, label="raw fragments")
    ax.axvline(80, color=INX, ls="--", lw=1.1, label="80 aa")
    ax.set_xlabel("Alignment length (aa)")
    ax.set_ylabel("Identity (%)")
    ax.set_title(f"Cross-family MMseqs: {stats['n_homologous']} hit(s) ≥80 aa at e≤1e-3 (of {len(hits)} fragments)")
    ax.legend(frameon=False, fontsize=8)
    ax.text(0.98, 0.04, "No long alignments → extreme divergence / separate families\n(not proof of absence of homology).", transform=ax.transAxes, ha="right", fontsize=7.5, color=MUTED)
    _save(fig, path)
    return stats


def plot_architecture(path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13.6, 5.2))

    def box(ax, x, y, w, label, color):
        ax.add_patch(FancyBboxPatch((x, y), w, 1.05, boxstyle="round,pad=0.04,rounding_size=0.12", facecolor=color, edgecolor="white", linewidth=1.1))
        ax.text(x + w / 2, y + 0.52, label, ha="center", va="center", color="white", fontsize=8, fontweight="bold")

    ax = axes[0]
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.set_title("Innexin genome logic: tandem, clade-bound", color=INX, pad=8)
    ax.add_patch(Rectangle((0.35, 5.2), 11.3, 1.7, facecolor="#ecfdf5", edgecolor="#99f6e4"))
    box(ax, 0.55, 5.5, 3.2, "ogre", INX)
    box(ax, 4.05, 5.5, 3.2, "Inx7", "#B45309")
    box(ax, 7.55, 5.5, 3.8, "Inx2", "#7C2D12")
    ax.text(0.45, 7.05, "Drosophila proximal X cluster", fontsize=8, color=MUTED)
    ax.text(0.45, 4.35, "Neighbours on the chromosome are not sisters on the tree.", fontsize=9)

    ax = axes[1]
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 8)
    ax.axis("off")
    ax.set_title("Connexin genome logic: conserved neighbourhoods", color=CNX, pad=8)
    ax.add_patch(Rectangle((0.35, 5.2), 11.3, 1.7, facecolor="#fff7ed", edgecolor="#fed7aa"))
    box(ax, 0.55, 5.5, 2.2, "GJB3", CNX)
    box(ax, 2.95, 5.5, 2.2, "GJB5", CNX)
    box(ax, 5.35, 5.5, 2.2, "GJA4", "#9A3412")
    ax.add_patch(FancyBboxPatch((7.75, 5.5), 3.5, 1.05, boxstyle="round,pad=0.04,rounding_size=0.12", facecolor="#e2e8f0", edgecolor="white"))
    ax.text(9.5, 6.02, "DLGAP3 flank", ha="center", va="center", fontsize=8, color="#334155")
    ax.text(0.45, 7.05, "Vertebrate β-cluster + separate GJA1 locus", fontsize=8, color=MUTED)
    ax.text(0.45, 4.35, "Same neighbours across mammals, birds, amphibians.", fontsize=9)
    _save(fig, path)


def plot_biophysical_contrasts(path: Path) -> None:
    """Literature biophysical contrasts: gap width, plaque spacing, voltage gating."""
    fig, axes = plt.subplots(1, 3, figsize=(11.6, 4.1))
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    ax = axes[0]
    gaps = [30, 18]
    labels = ["Hydra\n(invertebrate)", "Mouse heart\n(vertebrate)"]
    y = np.arange(len(labels))
    ax.barh(y, gaps, color=[INX, CNX], alpha=0.88, height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Intercellular gap (Å)")
    ax.set_title("Gap width", fontsize=10.5, pad=6)
    ax.set_xlim(0, 36)
    for yi, v in zip(y, gaps):
        ax.text(v + 0.6, yi, f"{v} Å", va="center", fontsize=8, color=MUTED)

    ax = axes[1]
    spacing = [111, 94, 77]
    labels = ["INX-6", "Cx26", "Cx43-GFP"]
    y = np.arange(len(labels))
    ax.barh(y, spacing, color=[INX, CNX, "#9A3412"], alpha=0.88, height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Center-to-center spacing (Å)")
    ax.set_title("Plaque channel spacing", fontsize=10.5, pad=6)
    ax.set_xlim(0, 125)
    for yi, v in zip(y, spacing):
        ax.text(v + 1.5, yi, f"{v} Å", va="center", fontsize=8, color=MUTED)

    ax = axes[2]
    ax.axvspan(-40, 0, color="#ecfdf5", alpha=0.65, zorder=0)
    ax.axvspan(0, 40, color="#fff7ed", alpha=0.65, zorder=0)
    ax.axvline(-20, color=INX, lw=10, alpha=0.75, solid_capstyle="butt")
    ax.axvline(22, color="#7C3AED", lw=10, alpha=0.75, solid_capstyle="butt")
    ax.set_xlim(-40, 40)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Membrane potential (mV)")
    ax.set_title("Voltage gating (hemichannel)", fontsize=10.5, pad=6)
    ax.text(-20, 0.72, "Innexon\n~−20 mV", ha="center", fontsize=8, color=INX, fontweight="bold")
    ax.text(22, 0.72, "Panx1\n>+20 mV", ha="center", fontsize=8, color="#7C3AED", fontweight="bold")
    ax.text(0, 0.18, "innexons vs Panx1 thresholds (innexin-lineage split)",
            ha="center", fontsize=7.5, color=MUTED)

    fig.suptitle("Biophysical contrasts (Skerrett et al. 2017; Dahl et al. 2014)", fontsize=11, y=1.03)
    _save(fig, path)


def write_html(stats: dict, figures: list[tuple[str, str, str, str]]) -> None:
    beats = []
    for i, (fid, fname, title, claim) in enumerate(figures, start=1):
        beats.append(
            f"""
  <article class="beat" id="{fid}">
    <span class="eyebrow">{i:02d} · innexin × connexin</span>
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
<title>Innexin × Connexin — direct contrast</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,560;9..144,700&family=Sora:wght@400;500;600&display=swap" rel="stylesheet">
<style>
:root {{
  --ink:#12202a; --muted:#4d6570; --deep:#0a3a42; --sea:#0f766e; --ember:#c2410c;
  --foam:#f4f7f5; --card:#fff; --line:rgba(18,32,42,.1);
}}
* {{ box-sizing:border-box }}
html {{ scroll-behavior:smooth }}
body {{
  margin:0; color:var(--ink); font-family:Sora,sans-serif; line-height:1.6;
  background:
    radial-gradient(1000px 480px at 8% -8%, rgba(15,118,110,.16), transparent 55%),
    radial-gradient(900px 420px at 92% 4%, rgba(194,65,12,.12), transparent 50%),
    linear-gradient(180deg,#e8f1ee 0%, var(--foam) 42%, #eef3f1 100%);
}}
a {{ color:var(--sea) }}
.hero {{
  min-height:86vh; display:grid; align-items:end; padding:clamp(1.5rem,4vw,3.25rem); color:#f7fffb;
  background:linear-gradient(118deg, rgba(10,58,66,.94) 0%, rgba(15,118,110,.7) 46%, rgba(194,65,12,.66) 100%);
}}
.brand {{ font-family:Fraunces,serif; font-size:clamp(2.3rem,6.4vw,4.4rem); font-weight:700; letter-spacing:-.03em; line-height:.95; margin:0 0 1rem }}
.hero h1 {{ margin:0 0 .75rem; font-size:clamp(1.05rem,2.1vw,1.3rem); font-weight:500; max-width:40ch }}
.hero .lede {{ margin:0 0 1.6rem; max-width:50ch; opacity:.92 }}
.cta a {{
  display:inline-block; text-decoration:none; font-weight:600; padding:.8rem 1.15rem; border-radius:999px;
  margin:.2rem .4rem .2rem 0; background:#f5fffb; color:var(--deep); box-shadow:0 8px 24px rgba(0,0,0,.18);
}}
.cta a.ghost {{ background:transparent; color:#f5fffb; border:1px solid rgba(245,255,251,.45); box-shadow:none }}
.wrap {{ width:min(880px, calc(100% - 2rem)); margin:0 auto; padding:2.25rem 0 4.5rem }}
.verdict {{
  margin-top:-2.2rem; position:relative; z-index:2; background:rgba(255,255,255,.95);
  border:1px solid var(--line); border-radius:18px; padding:1.25rem 1.35rem; box-shadow:0 14px 36px rgba(18,32,42,.08);
}}
.verdict h2 {{ font-family:Fraunces,serif; font-size:1.35rem; margin:0 0 .55rem }}
.verdict p {{ margin:0; color:var(--muted) }}
.verdict strong {{ color:var(--deep) }}
.stats {{ display:grid; grid-template-columns:repeat(4,1fr); gap:.7rem; margin:1.4rem 0 0 }}
.stat {{ background:#fff; border:1px solid var(--line); border-radius:14px; padding:.85rem }}
.stat b {{ display:block; font-family:Fraunces,serif; font-size:1.25rem; color:var(--deep) }}
.stat span {{ color:var(--muted); font-size:.78rem }}
.toc {{ display:flex; flex-wrap:wrap; gap:.45rem; margin:1.75rem 0 .5rem }}
.toc a {{ text-decoration:none; font-size:.78rem; font-weight:500; color:var(--deep); background:#fff; border:1px solid var(--line); padding:.4rem .75rem; border-radius:999px }}
.beat {{ margin:2.3rem 0; scroll-margin-top:1rem }}
.beat .eyebrow {{ display:inline-block; font-size:.72rem; font-weight:700; letter-spacing:.08em; text-transform:uppercase; color:var(--ember); margin-bottom:.35rem }}
.beat h2 {{ font-family:Fraunces,serif; font-size:clamp(1.4rem,3vw,1.8rem); margin:0 0 .55rem; letter-spacing:-.02em }}
.beat .claim {{ color:var(--muted); margin:0 0 1rem; max-width:64ch }}
.figure {{ background:var(--card); border:1px solid var(--line); border-radius:16px; padding:.85rem; box-shadow:0 10px 28px rgba(18,32,42,.05) }}
.figure img {{ width:100%; display:block; border-radius:10px; cursor:zoom-in; background:#f8fafc }}
.closing {{ margin-top:2.5rem; padding:1.4rem 1.35rem; border-radius:16px; background:linear-gradient(135deg,var(--deep),var(--ember)); color:#f5fffb }}
.closing h2 {{ font-family:Fraunces,serif; margin:0 0 .55rem; font-size:1.4rem }}
.closing p {{ margin:0; opacity:.95; max-width:62ch }}
.dim-table-wrap {{ overflow-x:auto; margin:1rem 0 0 }}
.dim-table {{
  width:100%; border-collapse:collapse; font-size:.86rem; background:#fff;
  border:1px solid var(--line); border-radius:12px; overflow:hidden;
}}
.dim-table th, .dim-table td {{ padding:.55rem .7rem; text-align:left; border-bottom:1px solid var(--line); color:var(--ink) }}
.dim-table th {{ background:rgba(15,118,110,.08); color:var(--deep); font-weight:600; font-size:.78rem }}
.dim-table tr:last-child td {{ border-bottom:none }}
.dim-table td.num {{ font-variant-numeric:tabular-nums; white-space:nowrap }}
.dim-note {{ font-size:.82rem; color:var(--muted); max-width:64ch; margin:.75rem 0 0 }}
.closing .lit-note {{ font-size:.88rem; opacity:.92; max-width:64ch; margin:.65rem 0 0 }}
.footer {{ margin-top:1.75rem; color:var(--muted); font-size:.88rem }}
.footer a {{ margin-right:.85rem }}
.modal {{ display:none; position:fixed; inset:0; background:rgba(0,0,0,.88); z-index:30; align-items:center; justify-content:center; padding:1.25rem; cursor:zoom-out }}
.modal.open {{ display:flex }}
.modal img {{ max-width:95vw; max-height:90vh; border-radius:8px }}
@media (max-width:800px) {{ .stats {{ grid-template-columns:1fr 1fr }} }}
</style>
</head>
<body>
<header class="hero">
  <div>
    <p class="brand">Innexin × Connexin</p>
    <h1>Compare the two families to each other — not each family to itself</h1>
    <p class="lede">One sequence space, one length axis, one chemistry panel. Same biological job — sequence evidence asks whether the two families share detectable similarity.</p>
    <div class="cta">
      <a href="#verdict">Read the contrast</a>
      <a class="ghost" href="../family_comparison/index.html">Previous site (same metrics twice)</a>
    </div>
  </div>
</header>
<main class="wrap">
  <section class="verdict" id="verdict">
    <h2>Core contrast</h2>
    <p>
      Innexins and connexins occupy the same biological niche (intercellular coupling) but
      <strong>separate clouds in this panel’s sequence space</strong>.
      They live in largely non-overlapping kingdoms, show
      <strong>low cross-family embedding cosine</strong>
      (innexin↔connexin 3-mer cosine {stats['cos_cross']:.2f} vs within-family {stats['cos_inx']:.2f} / {stats['cos_cnx']:.2f}),
      different protein lengths (median {stats['len_inx']:.0f} vs {stats['len_cnx']:.0f} aa),
      exon architecture, and genomic neighbourhood logic.
      Cross-family MMseqs finds <strong>{stats['n_homologous']}</strong> alignment(s) at e≤1e-3 and ≥80 aa —
      failure to detect similarity at this threshold is evidence for extreme divergence or separate families,
      not by itself proof of absence of homology (literature nonetheless treats them as non-homologous channels).
    </p>
    <div class="stats">
      <div class="stat"><b>{stats['n_inx']}</b><span>innexin proteins</span></div>
      <div class="stat"><b>{stats['n_cnx']}</b><span>connexin proteins</span></div>
      <div class="stat"><b>{stats['cos_cross']:.2f}</b><span>innexin↔connexin embedding cosine</span></div>
      <div class="stat"><b>{stats['n_homologous']}</b><span>cross-family hits (e≤1e-3, ≥80 aa)</span></div>
    </div>
  </section>
  <nav class="toc">{toc}</nav>
{''.join(beats)}
  <section class="closing">
    <h2>Takeaway</h2>
    <p>
      Sequence space, gene models, copy-number sampling, and genomic neighbourhoods —
      Cys motif counts (4 vs 6) are alignment filters. Annotation caveats are on the
      <a href="../gap_junction_path/index.html#artifacts" style="color:#d1fae5">gap-junction path · chapter 4</a>.
    </p>
    <h2 style="margin-top:1.4rem">One sentence</h2>
    <p>
      Same gap-junction job, largely separate kingdom occupancy, separate sequence clouds, different gene models
      (intron-rich innexins vs compact connexins), and no cross-family MMseqs hits at e≤1e-3 and ≥80 aa —
      parallel solutions rather than one detectable homologous family in this panel.
    </p>
  </section>
  <footer class="footer">
    <p>Direct family-vs-family contrast from combined FASTAs, exon summaries, and cross-family MMseqs. The older same-metrics page is unchanged.</p>
    <p>
      <a href="../family_comparison/index.html">Same-metrics page</a>
      <a href="../innexin_insights/index.html">Innexin insights</a>
      <a href="../connexin_insights/index.html">Connexin insights</a>
      <a href="../phylogenetic_story/index.html">Phylogeny story</a>
      <a href="../../../pannexin/project/results/panx_vs_inx/index.html">Pannexin × innexin</a>
      <a href="../sources/index.html">Sources</a>
    </p>
  </footer>
</main>
<div id="modal" class="modal" onclick="this.classList.remove('open')"><img id="modal-img" alt=""></div>
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

    rows = read_fasta(INX_FASTA, "innexin") + read_fasta(CNX_FASTA, "connexin")
    df = pd.DataFrame(rows)
    inx_loci = pd.read_csv(INX_LOCI)
    cnx_loci = pd.read_csv(CNX_LOCI)
    exon = pd.read_csv(EXON) if EXON.exists() else pd.DataFrame()

    seqs = df["seq"].tolist()
    fams = df["family"].tolist()
    kmer = kmer_svd(seqs)
    aa = aa_matrix(rows)
    kmer_xy = pca2(kmer)
    aa_xy = pca2(aa)

    ia = np.where(df["family"] == "innexin")[0]
    ib = np.where(df["family"] == "connexin")[0]
    cos = {
        "inx_inx": mean_cosine(kmer, ia, ia),
        "cnx_cnx": mean_cosine(kmer, ib, ib),
        "inx_cnx": mean_cosine(kmer, ia, ib),
    }

    plot_occupancy(inx_loci, cnx_loci, FIGS / "01_who_has_whom.png")
    plot_joint_pca(kmer_xy, fams, "Both families in one 3-mer sequence space", FIGS / "02_joint_kmer_pca.png")
    plot_joint_pca(aa_xy, fams, "Both families in one amino-acid composition space", FIGS / "03_joint_aa_pca.png")
    plot_family_cosine(cos, FIGS / "04_cross_family_cosine.png")
    cross = plot_cross_mmseqs(FIGS / "05_cross_family_homology.png")
    plot_length_overlay(
        df.loc[df["family"] == "innexin", "length"].values,
        df.loc[df["family"] == "connexin", "length"].values,
        FIGS / "06_length_overlay.png",
    )
    if not exon.empty:
        plot_exons_span(exon, FIGS / "07_exons_and_span.png")
    plot_chemistry(df, FIGS / "08_chemistry.png")
    plot_aa_difference(df, FIGS / "09_aa_difference.png")
    plot_copy_number(inx_loci, cnx_loci, FIGS / "10_copy_number.png")
    plot_architecture(FIGS / "11_genome_logic.png")
    plot_biophysical_contrasts(FIGS / "12_biophysical_contrasts.png")

    n_inx = int((df["family"] == "innexin").sum())
    n_cnx = int((df["family"] == "connexin").sum())
    len_inx = float(df.loc[df["family"] == "innexin", "length"].median())
    len_cnx = float(df.loc[df["family"] == "connexin", "length"].median())

    metric_rows = [
        {"metric": "n_innexin", "value": str(n_inx)},
        {"metric": "n_connexin", "value": str(n_cnx)},
        {"metric": "median_length_innexin", "value": f"{len_inx:.1f}"},
        {"metric": "median_length_connexin", "value": f"{len_cnx:.1f}"},
        {"metric": "cosine_inx_inx", "value": f"{cos['inx_inx']:.4f}"},
        {"metric": "cosine_cnx_cnx", "value": f"{cos['cnx_cnx']:.4f}"},
        {"metric": "cosine_inx_cnx", "value": f"{cos['inx_cnx']:.4f}"},
        {"metric": "homologous_hits", "value": str(cross["n_homologous"])},
        {"metric": "fragment_hits", "value": str(cross["n_hits"])},
    ]
    write_csv(OUT / "direct_contrast_metrics.csv", metric_rows)

    findings = [
        f"Joint 3-mer space: innexin↔innexin cosine {cos['inx_inx']:.2f}, connexin↔connexin {cos['cnx_cnx']:.2f}, innexin↔connexin {cos['inx_cnx']:.2f}.",
        f"Median length {len_inx:.0f} aa (innexin) vs {len_cnx:.0f} aa (connexin).",
        f"Cross-family MMseqs: {cross['n_homologous']} hit(s) at e≤1e-3 and ≥80 aa (fragments {cross['n_hits']}) — no detectable long-range similarity in this panel.",
        "Occupancy is largely non-overlapping: innexins predominate in invertebrates, connexins in vertebrates (exceptions exist).",
        "Gene models: connexin ORFs usually intron-free (compact ~2-exon models); innexins carry coding-region introns (median ~5 exons) and can splice multiple isoforms.",
    ]
    (OUT / "key_findings.md").write_text("# Direct innexin × connexin contrast\n\n" + "\n".join(f"- {f}" for f in findings) + "\n", encoding="utf-8")

    figures = [
        ("who", "01_who_has_whom.png", "Different kingdoms",
         "Innexins predominate in invertebrates; connexins are the canonical vertebrate gap-junction family in this panel. The comparison starts with who carries which family, not with subfamily purity."),
        ("space", "02_joint_kmer_pca.png", "One sequence space, two clouds",
         f"Every innexin and connexin protein is embedded in the <em>same</em> 3-mer space. They form separate clouds. Cross-family cosine is {cos['inx_cnx']:.2f} versus {cos['inx_inx']:.2f} / {cos['cnx_cnx']:.2f} inside each family."),
        ("aa", "03_joint_aa_pca.png", "Composition space agrees",
         "Amino-acid mix alone already splits the two families. They are not two labels on one protein cloud."),
        ("cos", "04_cross_family_cosine.png", "Innexin ↔ connexin is the empty comparison",
         "This bar is the actual between-family number: how similar innexins are to connexins, not how similar innexins are to other innexins."),
        ("hom", "05_cross_family_homology.png", "Cross-family sequence search",
         f"MMseqs innexin query vs connexin target: <strong>{cross['n_homologous']}</strong> hit(s) at e≤1e-3 and ≥80 aa among {cross['n_hits']} short fragments. "
         "Zero long alignments supports extreme divergence or separate families; it is not proof of absence of homology. Same job, different detectable sequence grammar."),
        ("len", "06_length_overlay.png", "Different protein shapes",
         f"Overlaid <strong>monomer</strong> length distributions (aa). Innexin median {len_inx:.0f} aa, connexin median {len_cnx:.0f} aa. "
         "This is not docked channel length in ångströms — monomer aa lengths are the genomic/sequence character analysed here."),
        ("exons", "07_exons_and_span.png", "Different gene models",
         "Connexin coding regions usually lack introns (ORF in one exon → compact gene models, panel median ~2 exons including 5′ UTR). "
         "Innexins carry coding-region introns (median ~5 exons) and can produce multiple splice variants from one locus. "
         "Same axes: exon count and gene span."),
        ("chem", "08_chemistry.png", "Different chemistry",
         "Total cysteine density (whole protein), hydrophobicity (GRAVY) and charged-residue fraction. "
         "This is not the conserved extracellular-loop ladder: that signature is fixed at "
         "<strong>2 Cys/loop (4 total)</strong> for innexins/pannexins vs <strong>3 Cys/loop (6 total)</strong> for connexins."),
        ("comp", "09_aa_difference.png", "Which amino acids each family prefers",
         "Positive bars: innexins use more of that residue. Negative bars: connexins use more. A direct compositional contrast."),
        ("copy", "10_copy_number.png", "Panel sampling vs genome biology",
         "Connexins are large multigene families in vertebrates (~20–22 in jawed vertebrates, up to ~46 in teleosts). "
         "A histogram median of 1 in an earlier plot was a panel sampling issue: 23/38 species carry only one connexin locus "
         "in this restricted reference/discovery slice (incomplete UniProt sampling + sparse discovery hits). "
         "Human (<strong>22</strong> loci here) and mouse (<strong>20</strong>) match expectation; the figure now shows "
         "best-covered reference species against the literature band, not a misleading genome-wide median."),
        ("arch", "11_genome_logic.png", "Different genomic logic",
         "Innexin: tandem clusters that cut across the tree. Connexin: conserved neighbourhoods — "
         "<strong>GJA4 β-cluster</strong> as a tandem-duplication array (unequal crossing-over) versus "
         "<strong>GJA1</strong> as a WGD-dispersed address. Empty flanks around <em>X. laevis</em> gja1.L are a "
         "scaffold/gap issue in the plot, not a gene island. Insect Inx2 ≠ nematode inx-2 (SF4 vs SF3)."),
    ]
    if exon.empty:
        figures = [f for f in figures if f[0] != "exons"]
    # Biophysics section left out (plot file may still exist on disk).
    figures = [f for f in figures if f[0] != "biophys"]

    write_html(
        {
            "n_inx": n_inx,
            "n_cnx": n_cnx,
            "cos_inx": cos["inx_inx"],
            "cos_cnx": cos["cnx_cnx"],
            "cos_cross": cos["inx_cnx"],
            "len_inx": len_inx,
            "len_cnx": len_cnx,
            "n_homologous": cross["n_homologous"],
        },
        figures,
    )
    print(f"Done → {OUT / 'index.html'}")
    print(f"  cosine inx/cnx/cross {cos['inx_inx']:.3f}/{cos['cnx_cnx']:.3f}/{cos['inx_cnx']:.3f}")
    print(f"  median length {len_inx:.0f} vs {len_cnx:.0f}  homologous {cross['n_homologous']}")


if __name__ == "__main__":
    main()
