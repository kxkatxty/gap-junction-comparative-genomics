#!/usr/bin/env python3
"""Guided phylogenetic story path for gap-junction proteins.

New hub page that explains how innexins/connexins work, then sends the reader
into the existing summed pages (insights, clade comparisons, inx_vs_cnx).

Does not rebuild those sites. Outputs → project/results/gap_junction_path/
"""

from __future__ import annotations

from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

from pipeline.common import RESULTS_DIR

OUT = RESULTS_DIR / "gap_junction_path"
FIGS = OUT / "figures"
INX = "#0F766E"
CNX = "#C2410C"
INK = "#12202a"
MUTED = "#4d6570"
MEM = "#1e3a4c"

plt.rcParams.update({"figure.dpi": 150, "savefig.dpi": 200, "font.size": 10})


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_two_cells(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(12.4, 5.5))
    ax.set_xlim(0, 14)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("Gap junctions couple neighbouring cells", loc="left", color=INK, fontsize=13, pad=8)

    ax.add_patch(FancyBboxPatch((0.3, 0.7), 4.55, 4.6, boxstyle="round,pad=0.08,rounding_size=0.45", facecolor="#dbece8", edgecolor=INX, linewidth=1.6))
    ax.add_patch(FancyBboxPatch((9.15, 0.7), 4.55, 4.6, boxstyle="round,pad=0.08,rounding_size=0.45", facecolor="#fde8dc", edgecolor=CNX, linewidth=1.6))
    ax.text(2.57, 4.85, "Cell A", ha="center", fontsize=11, color=INX, fontweight="bold")
    ax.text(11.42, 4.85, "Cell B", ha="center", fontsize=11, color=CNX, fontweight="bold")

    # membranes facing the wider extracellular gap
    ax.add_patch(Rectangle((4.63, 0.85), 0.22, 4.3, facecolor=MEM, alpha=0.85, zorder=2))
    ax.add_patch(Rectangle((9.15, 0.85), 0.22, 4.3, facecolor=MEM, alpha=0.85, zorder=2))
    ax.text(7.0, 0.28, "extracellular gap  (~2–4 nm)", ha="center", fontsize=8, color=MUTED)

    # channel spanning the gap
    ax.add_patch(FancyBboxPatch((4.35, 2.2), 5.3, 1.65, boxstyle="round,pad=0.02,rounding_size=0.4", facecolor="#f8fafc", edgecolor=MEM, linewidth=1.4, zorder=3))
    ax.add_patch(Rectangle((4.85, 2.45), 4.3, 1.15, facecolor="#94a3b8", alpha=0.28, zorder=4))
    ax.text(7.0, 3.02, "gap junction proteins", ha="center", va="center", fontsize=10, color=INK, fontweight="bold", zorder=5)

    for x, col in ((1.35, INX), (2.15, "#0e7490"), (2.95, "#115e59")):
        ax.add_patch(Circle((x, 3.0), 0.13, color=col, zorder=5, alpha=0.85))
    for x, col in ((11.05, CNX), (11.85, "#9a3412"), (12.65, "#c2410c")):
        ax.add_patch(Circle((x, 3.0), 0.13, color=col, zorder=5, alpha=0.85))
    ax.annotate("", xy=(4.45, 3.02), xytext=(3.25, 3.02), arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.4))
    ax.annotate("", xy=(10.85, 3.02), xytext=(9.55, 3.02), arrowprops=dict(arrowstyle="->", color=MUTED, lw=1.4))
    ax.text(2.55, 1.35, "ions · IP3 · cAMP\nsmall metabolites", ha="center", fontsize=8, color=MUTED)
    ax.text(11.4, 1.35, "electrical + metabolic\ncoupling", ha="center", fontsize=8, color=MUTED)
    _save(fig, path)


def plot_topology(path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(11.4, 5.0))
    for ax, title, color, n_cys in (
        (axes[0], "Innexin monomer  (4 TM)", INX, 2),
        (axes[1], "Connexin monomer  (4 TM)", CNX, 3),
    ):
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 8)
        ax.axis("off")
        ax.set_title(title, color=color, pad=8)
        ax.axhspan(3.15, 4.85, color="#cbd5e1", alpha=0.55)
        ax.text(0.35, 5.05, "extracellular", fontsize=8, color=MUTED)
        ax.text(0.35, 2.55, "cytoplasm", fontsize=8, color=MUTED)
        xs = [1.6, 3.4, 5.2, 7.0]
        for i, x in enumerate(xs, start=1):
            ax.add_patch(FancyBboxPatch((x, 3.2), 1.15, 1.6, boxstyle="round,pad=0.02,rounding_size=0.12", facecolor=color, edgecolor="white", linewidth=1.2))
            ax.text(x + 0.57, 4.0, f"TM{i}", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
        # EL1: TM1–TM2 extracellular
        ax.plot([2.17, 2.17], [4.8, 6.35], color=color, lw=2)
        ax.plot([2.17, 3.97], [6.35, 6.35], color=color, lw=2)
        ax.plot([3.97, 3.97], [6.35, 4.8], color=color, lw=2)
        ax.scatter(np.linspace(2.45, 3.7, n_cys), [6.35] * n_cys, s=28, c="#fbbf24", zorder=5, edgecolors="white")
        # cytoplasmic loop TM2–TM3
        ax.plot([3.97, 3.97], [3.2, 1.65], color=color, lw=2)
        ax.plot([3.97, 5.77], [1.65, 1.65], color=color, lw=2)
        ax.plot([5.77, 5.77], [1.65, 3.2], color=color, lw=2)
        # EL2: TM3–TM4 extracellular
        ax.plot([5.77, 5.77], [4.8, 6.35], color=color, lw=2)
        ax.plot([5.77, 7.57], [6.35, 6.35], color=color, lw=2)
        ax.plot([7.57, 7.57], [6.35, 4.8], color=color, lw=2)
        ax.scatter(np.linspace(6.05, 7.3, n_cys), [6.35] * n_cys, s=28, c="#fbbf24", zorder=5, edgecolors="white")
        ax.plot([1.6, 0.7], [3.2, 2.3], color=color, lw=2)
        ax.plot([8.15, 9.2], [3.2, 2.3], color=color, lw=2)
        ax.text(0.55, 1.95, "N", color=color, fontsize=10, fontweight="bold")
        ax.text(9.05, 1.95, "C", color=color, fontsize=10, fontweight="bold")
        ax.text(4.9, 6.85, f"conserved Cys  ({n_cys} per loop · {n_cys * 2} total)", ha="center", fontsize=7.5, color="#b45309")
        ax.text(5.0, 0.45, "Same 4-TM topology. Different Cys ladder & sequence family.", ha="center", fontsize=8, color=MUTED)
    _save(fig, path)


def plot_hexamers(path: Path) -> None:
    """Innexon = octamer (8); connexon = hexamer (6). Kept filename for existing HTML links."""
    fig, ax = plt.subplots(figsize=(10.8, 5.8))
    ax.set_xlim(-1.2, 11.2)
    ax.set_ylim(-0.5, 6.3)
    ax.axis("off")
    ax.set_title(
        "Different stoichiometry, same docking idea",
        loc="left",
        color=INK,
        fontsize=13,
        pad=6,
    )

    def ring(cx, cy, n, color, label, radius=1.05, bead=0.38):
        for k in range(n):
            ang = np.pi / 2 + k * (2 * np.pi / n)
            x = cx + radius * np.cos(ang)
            y = cy + radius * np.sin(ang)
            ax.add_patch(Circle((x, y), bead, facecolor=color, edgecolor="white", lw=1.3, alpha=0.92))
        ax.add_patch(Circle((cx, cy), 0.42, facecolor="#f8fafc", edgecolor=color, lw=1.3))
        ax.text(cx, cy, "pore", ha="center", va="center", fontsize=7.5, color=MUTED)
        ax.text(cx, cy - 1.85, label, ha="center", fontsize=9, color=color, fontweight="bold")

    # Cryo-EM: innexin hemichannels are octameric; connexin hemichannels are hexameric.
    ring(2.25, 3.55, 8, INX, "innexon  (innexin octamer · 8)", radius=1.12, bead=0.34)
    ring(8.75, 3.55, 6, CNX, "connexon  (connexin hexamer · 6)", radius=0.95, bead=0.40)
    ax.annotate("", xy=(6.85, 3.55), xytext=(4.2, 3.55), arrowprops=dict(arrowstyle="<->", color=MUTED, lw=1.6))
    ax.text(5.5, 4.35, "dock across\nthe gap", ha="center", fontsize=8, color=MUTED)
    ax.text(
        5.5,
        0.45,
        "Complete gap junction: innexin 8 + 8 = 16 subunits · connexin 6 + 6 = 12 subunits.",
        ha="center",
        fontsize=9,
        color=INK,
    )
    ax.text(
        5.5,
        0.05,
        "Innexins therefore have a larger outer diameter and a wider pore.",
        ha="center",
        fontsize=8,
        color=MUTED,
    )
    _save(fig, path)


def plot_kingdoms(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(10.6, 4.8))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6)
    ax.axis("off")
    ax.set_title("Same job, two evolutionary experiments", loc="left", color=INK, fontsize=13, pad=6)
    fork_y = 3.45
    box_top = 2.12
    ax.plot([6, 6], [5.3, fork_y], color="#94a3b8", lw=2, zorder=1)
    ax.plot([6, 2.85], [fork_y, box_top], color=INX, lw=2.2, zorder=1, solid_capstyle="round")
    ax.plot([6, 9.15], [fork_y, box_top], color=CNX, lw=2.2, zorder=1, solid_capstyle="round")
    ax.add_patch(Circle((6, 5.35), 0.22, color="#334155", zorder=2))
    ax.text(6, 5.75, "Metazoa", ha="center", fontsize=9, color=INK)
    ax.add_patch(FancyBboxPatch((0.55, 0.4), 4.6, 1.62, boxstyle="round,pad=0.0,rounding_size=0.16", facecolor="#ecfdf5", edgecolor=INX, lw=1.5, zorder=3))
    ax.add_patch(FancyBboxPatch((6.85, 0.4), 4.6, 1.62, boxstyle="round,pad=0.0,rounding_size=0.16", facecolor="#fff7ed", edgecolor=CNX, lw=1.5, zorder=3))
    ax.add_patch(Circle((2.85, box_top), 0.09, color=INX, zorder=4))
    ax.add_patch(Circle((9.15, box_top), 0.09, color=CNX, zorder=4))
    ax.text(2.85, 1.42, "Innexins", ha="center", fontsize=13, color=INX, fontweight="bold", zorder=4)
    ax.text(2.85, 0.82, "insects · nematodes · molluscs\nrotifers · other invertebrates", ha="center", fontsize=8, color=MUTED, zorder=4)
    ax.text(9.15, 1.42, "Connexins", ha="center", fontsize=13, color=CNX, fontweight="bold", zorder=4)
    ax.text(9.15, 0.82, "mammals · birds · amphibians\nfish · other vertebrates", ha="center", fontsize=8, color=MUTED, zorder=4)
    ax.text(6, 3.68, "separate sequence families (lit.)", ha="center", fontsize=7.5, color=MUTED, style="italic")
    _save(fig, path)


def plot_more_than_pipes(path: Path) -> None:
    """Beyond intercellular pipes: lipid seal, morphogenetic scaffold, non-junctional ATP."""
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.9))
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 8)
        ax.axis("off")

    # Panel 1: lipid-mediated N-terminal gating
    ax = axes[0]
    ax.set_title("N-terminal lipid gating", fontsize=11, color=INK, pad=6)
    ax.add_patch(Rectangle((0, 3.2), 10, 1.6, color="#dbe4ea"))
    ax.add_patch(Circle((5, 4.0), 1.15, facecolor=INX, edgecolor="#115e59", lw=1.4, alpha=0.9))
    ax.add_patch(Circle((5, 4.0), 0.55, facecolor="#f8fafc", edgecolor="#115e59", lw=1))
    # N-termini (flexible) + double-layer lipid seal
    for dx in (-0.95, -0.35, 0.35, 0.95):
        ax.plot([5 + dx, 5 + dx * 0.55], [5.05, 4.55], color="#115e59", lw=1.6, solid_capstyle="round", zorder=4)
    ax.add_patch(FancyBboxPatch((3.85, 3.55), 2.3, 0.35, boxstyle="round,pad=0.01,rounding_size=0.08",
                                facecolor="#fbbf24", edgecolor="#b45309", lw=1, zorder=5))
    ax.add_patch(FancyBboxPatch((3.85, 4.05), 2.3, 0.35, boxstyle="round,pad=0.01,rounding_size=0.08",
                                facecolor="#fde68a", edgecolor="#b45309", lw=1, zorder=5))
    ax.text(5, 4.0, "lipids", ha="center", va="center", fontsize=7.5, color="#7c2d12", fontweight="bold", zorder=6)
    ax.text(5, 1.9, "INX-6 in nanodisc (cryo-EM)\nN-term rearrangement → lipid seal", ha="center", fontsize=8, color=MUTED)
    ax.text(5, 0.55, "active gating · intact N-term", ha="center", fontsize=8.5, color=INX, fontweight="bold")

    # Panel 2: Inx3 morphogenetic hierarchy
    ax = axes[1]
    ax.set_title("Inx3 morphogenetic hierarchy", fontsize=11, color=INK, pad=6)
    ax.add_patch(Rectangle((0, 3.4), 10, 1.3, color="#dbe4ea"))
    ax.plot([1.5, 8.5], [4.05, 4.05], color="#94a3b8", lw=6, solid_capstyle="round", alpha=0.5)
    # Inx3 hub (upstream)
    ax.add_patch(FancyBboxPatch((3.85, 5.05), 2.3, 0.72, boxstyle="round,pad=0.02,rounding_size=0.12",
                                facecolor="#0e7490", edgecolor="white", lw=1.2))
    ax.text(5, 5.41, "Inx3", ha="center", va="center", color="white", fontsize=9, fontweight="bold")
    # Dependent partners at the membrane
    ax.add_patch(FancyBboxPatch((1.4, 3.55), 1.55, 0.95, boxstyle="round,pad=0.02,rounding_size=0.12",
                                facecolor=INX, edgecolor="white", lw=1))
    ax.add_patch(FancyBboxPatch((4.2, 3.55), 1.55, 0.95, boxstyle="round,pad=0.02,rounding_size=0.12",
                                facecolor=INX, edgecolor="white", lw=1))
    ax.add_patch(FancyBboxPatch((7.0, 3.55), 1.55, 0.95, boxstyle="round,pad=0.02,rounding_size=0.12",
                                facecolor="#c2410c", edgecolor="white", lw=1))
    ax.text(2.18, 4.02, "Inx1", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax.text(5.0, 4.02, "Inx2", ha="center", va="center", color="white", fontsize=8, fontweight="bold")
    ax.text(7.78, 4.02, "DE-cad", ha="center", va="center", color="white", fontsize=7.5, fontweight="bold")
    for tx in (2.18, 5.0, 7.78):
        ax.annotate("", xy=(tx, 4.55), xytext=(5, 5.05),
                    arrowprops=dict(arrowstyle="->", color="#0e7490", lw=1.2))
    ax.text(5, 1.9, "Dorsal closure: Inx3 upstream\nof Inx1/Inx2 · parallel to DE-cad", ha="center", fontsize=8, color=MUTED)
    ax.text(5, 0.55, "Df-Inx3 → PM loss · lysosomal degradation", ha="center", fontsize=7.8, color="#0e7490", fontweight="bold")

    # Panel 3: non-junctional ATP
    ax = axes[2]
    ax.set_title("Non-junctional signaling", fontsize=11, color=INK, pad=6)
    ax.add_patch(Rectangle((0, 2.8), 10, 1.5, color="#dbe4ea"))
    ax.add_patch(Circle((5, 3.55), 0.95, facecolor="#c2410c", edgecolor="#7c2d12", lw=1.4, alpha=0.92))
    ax.add_patch(Circle((5, 3.55), 0.32, facecolor="#f8fafc", edgecolor="#7c2d12", lw=1))
    ax.annotate("ATP", xy=(5, 4.6), xytext=(5, 6.2),
                arrowprops=dict(arrowstyle="->", color="#c2410c", lw=1.5),
                ha="center", fontsize=10, color="#c2410c", fontweight="bold")
    ax.text(8.2, 5.5, "go", ha="center", fontsize=8, color="#c2410c",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#fff7ed", edgecolor="#fdba74"))
    ax.text(8.2, 4.5, "NO\nwhere", ha="center", fontsize=7.5, color="#0369a1",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#eff6ff", edgecolor="#93c5fd"))
    ax.text(8.2, 3.2, "AA\nstop", ha="center", fontsize=7.5, color="#7c2d12",
            bbox=dict(boxstyle="round,pad=0.2", facecolor="#fef2f2", edgecolor="#fca5a5"))
    ax.text(5, 1.55, "Innexon / pannexon\nATP release (e.g. neuroglia)", ha="center", fontsize=8, color=MUTED)
    ax.text(5, 0.45, "channel, not only pipe", ha="center", fontsize=8.5, color="#c2410c", fontweight="bold")

    fig.suptitle("Gap-junction proteins do more than couple two cytoplasms", fontsize=12, color=INK, y=1.02)
    _save(fig, path)


def plot_biophysical_contrasts(path: Path) -> None:
    """Literature biophysical contrasts: gap width, plaque spacing, voltage gating."""
    fig, axes = plt.subplots(1, 3, figsize=(11.6, 4.1))
    for ax in axes:
        ax.spines["top"].set_visible(False)
        ax.spines["right"].set_visible(False)

    # Intercellular gap width
    ax = axes[0]
    gaps = [30, 18]
    labels = ["Hydra\n(invertebrate)", "Mouse heart\n(vertebrate)"]
    y = np.arange(len(labels))
    ax.barh(y, gaps, color=[INX, CNX], alpha=0.88, height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Intercellular gap (Å)")
    ax.set_title("Gap width", fontsize=10.5, color=INK, pad=6)
    ax.set_xlim(0, 36)
    for yi, v in zip(y, gaps):
        ax.text(v + 0.6, yi, f"{v} Å", va="center", fontsize=8, color=MUTED)

    # Plaque center-to-center spacing
    ax = axes[1]
    spacing = [111, 94, 77]
    labels = ["INX-6", "Cx26", "Cx43-GFP"]
    y = np.arange(len(labels))
    ax.barh(y, spacing, color=[INX, CNX, "#9A3412"], alpha=0.88, height=0.58)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8.5)
    ax.set_xlabel("Center-to-center spacing (Å)")
    ax.set_title("Plaque channel spacing", fontsize=10.5, color=INK, pad=6)
    ax.set_xlim(0, 125)
    for yi, v in zip(y, spacing):
        ax.text(v + 1.5, yi, f"{v} Å", va="center", fontsize=8, color=MUTED)

    # Voltage gating thresholds (hemichannel)
    ax = axes[2]
    ax.axvspan(-40, 0, color="#ecfdf5", alpha=0.65, zorder=0)
    ax.axvspan(0, 40, color="#fff7ed", alpha=0.65, zorder=0)
    ax.axvline(-20, color=INX, lw=10, alpha=0.75, solid_capstyle="butt")
    ax.axvline(22, color="#7C3AED", lw=10, alpha=0.75, solid_capstyle="butt")
    ax.set_xlim(-40, 40)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.set_xlabel("Membrane potential (mV)")
    ax.set_title("Voltage gating (hemichannel)", fontsize=10.5, color=INK, pad=6)
    ax.text(-20, 0.72, "Innexon\n~−20 mV", ha="center", fontsize=8, color=INX, fontweight="bold")
    ax.text(22, 0.72, "Panx1\n>+20 mV", ha="center", fontsize=8, color="#7C3AED", fontweight="bold")
    ax.text(0, 0.18, "mild depolarization opens innexons · Panx1 needs unphysiological positive Vm",
            ha="center", fontsize=7.5, color=MUTED)

    fig.suptitle("Biophysical contrasts at gap junctions (Skerrett et al. 2017; Dahl et al. 2014)", fontsize=11, color=INK, y=1.03)
    _save(fig, path)


def plot_chordate_bottleneck(path: Path) -> None:
    """Literature scenario: early chordate loss of innexin diversity (Welzel & Schuster 2022)."""
    fig, ax = plt.subplots(figsize=(11.0, 4.6))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 6.2)
    ax.axis("off")
    ax.set_title(
        "Early chordate bottleneck → connexins (literature sketch)",
        loc="left",
        color=INK,
        fontsize=13,
        pad=6,
    )

    # timeline
    ax.plot([0.8, 11.2], [3.35, 3.35], color="#94a3b8", lw=2.2, zorder=1)
    stages = [
        (2.0, INX, "Non-chordates", "many innexins\n± NGS · rich GJ diversity"),
        (6.0, "#b45309", "Lancelets\n(basal chordates)", "1 conserved innexin\nwith NGS · no connexins"),
        (10.0, CNX, "Jawed vertebrates", "connexins diversify\n+ pannexins (NGS)"),
    ]
    for x, color, title, note in stages:
        ax.add_patch(Circle((x, 3.35), 0.18, facecolor=color, edgecolor="white", lw=1.5, zorder=3))
        ax.add_patch(
            FancyBboxPatch(
                (x - 1.55, 0.45),
                3.1,
                2.15,
                boxstyle="round,pad=0.02,rounding_size=0.14",
                facecolor="white",
                edgecolor=color,
                lw=1.6,
                zorder=2,
            )
        )
        ax.text(x, 2.05, title, ha="center", va="center", fontsize=10, fontweight="bold", color=color, zorder=3)
        ax.text(x, 1.15, note, ha="center", va="center", fontsize=8, color=MUTED, zorder=3)

    ax.annotate(
        "",
        xy=(4.35, 3.35),
        xytext=(3.65, 3.35),
        arrowprops=dict(arrowstyle="->", color="#94a3b8", lw=1.8),
    )
    ax.annotate(
        "",
        xy=(8.35, 3.35),
        xytext=(7.65, 3.35),
        arrowprops=dict(arrowstyle="->", color="#94a3b8", lw=1.8),
    )
    ax.text(4.0, 3.85, "diversity collapse", ha="center", fontsize=8, color="#b45309", style="italic")
    ax.text(8.0, 3.85, "de novo Cx family", ha="center", fontsize=8, color=CNX, style="italic")

    ax.text(
        6.0,
        5.35,
        "Literature hypothesis (Welzel & Schuster 2022, eLife): loss of innexin diversity at the chordate origin\n"
        "— especially one glycosylated innexin in lancelets — as pressure for connexins to restore intercellular pathways.",
        ha="center",
        fontsize=8.5,
        color=MUTED,
    )
    ax.text(
        6.0,
        4.55,
        "Not measured in this panel (no Amphioxus search here) — evolutionary context for the binary story.",
        ha="center",
        fontsize=8,
        color="#64748b",
    )
    _save(fig, path)


def write_html() -> None:
    html = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>A phylogenetic story of gap junctions</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,560;9..144,700&family=Sora:wght@400;500;600&display=swap" rel="stylesheet">
<style>
:root {
  --ink:#12202a; --muted:#4d6570; --deep:#0a3a42; --sea:#0f766e; --ember:#c2410c;
  --foam:#f4f7f5; --card:#fff; --line:rgba(18,32,42,.1);
}
* { box-sizing:border-box }
html { scroll-behavior:smooth }
body {
  margin:0; color:var(--ink); font-family:Sora,sans-serif; line-height:1.65;
  background:
    radial-gradient(1000px 480px at 8% -8%, rgba(15,118,110,.16), transparent 55%),
    radial-gradient(900px 420px at 92% 4%, rgba(194,65,12,.12), transparent 50%),
    linear-gradient(180deg,#e8f1ee 0%, var(--foam) 40%, #eef3f1 100%);
}
a { color:var(--sea) }
.hero {
  min-height:88vh; display:grid; align-items:end; padding:clamp(1.4rem,4vw,3.2rem); color:#f7fffb;
  background:linear-gradient(118deg, rgba(10,58,66,.94) 0%, rgba(15,118,110,.7) 48%, rgba(194,65,12,.62) 100%);
}
.brand { font-family:Fraunces,serif; font-size:clamp(2.3rem,6.5vw,4.5rem); font-weight:700; letter-spacing:-.03em; line-height:.95; margin:0 0 1rem }
.hero h1 { margin:0 0 .7rem; font-size:clamp(1.05rem,2.1vw,1.28rem); font-weight:500; max-width:38ch }
.hero .lede { margin:0 0 1.5rem; max-width:50ch; opacity:.93 }
.cta a {
  display:inline-block; text-decoration:none; font-weight:600; padding:.8rem 1.15rem; border-radius:999px;
  margin:.2rem .35rem .2rem 0; background:#f5fffb; color:var(--deep);
}
.cta a.ghost { background:transparent; color:#f5fffb; border:1px solid rgba(245,255,251,.45) }
.path {
  position:sticky; top:0; z-index:8; background:rgba(244,247,245,.92); backdrop-filter:blur(10px);
  border-bottom:1px solid var(--line); padding:.65rem 1rem;
}
.path-inner { width:min(920px,100%); margin:0 auto; display:flex; flex-wrap:wrap; gap:.4rem }
.path a {
  text-decoration:none; font-size:.74rem; font-weight:600; color:var(--deep); background:#fff;
  border:1px solid var(--line); padding:.35rem .65rem; border-radius:999px;
}
.path a:hover { border-color:var(--sea); color:var(--sea) }
.wrap { width:min(860px, calc(100% - 2rem)); margin:0 auto; padding:2rem 0 4.5rem }
.chapter { margin:2.4rem 0; scroll-margin-top:3.4rem }
.eyebrow { display:inline-block; font-size:.72rem; font-weight:700; letter-spacing:.08em; text-transform:uppercase; color:var(--ember); margin-bottom:.35rem }
.chapter h2 { font-family:Fraunces,serif; font-size:clamp(1.45rem,3vw,1.95rem); margin:0 0 .6rem; letter-spacing:-.02em }
.chapter p { color:var(--muted); margin:0 0 1rem; max-width:64ch }
.figure { background:var(--card); border:1px solid var(--line); border-radius:16px; padding:.85rem; margin:0 0 1rem; box-shadow:0 10px 28px rgba(18,32,42,.05) }
.figure img { width:100%; display:block; border-radius:10px; cursor:zoom-in; background:#f8fafc }
.cap { margin:.65rem .1rem 0; font-size:.88rem; color:var(--muted) }
.take { margin:.2rem 0 0; padding:.75rem .9rem; border-left:3px solid var(--sea); background:rgba(15,118,110,.07); border-radius:0 12px 12px 0; font-size:.92rem; color:var(--deep) }
.cards { display:grid; gap:.8rem; margin:1rem 0 0 }
.card {
  display:block; text-decoration:none; color:inherit; background:#fff; border:1px solid var(--line);
  border-radius:16px; padding:1rem 1.1rem; box-shadow:0 8px 22px rgba(18,32,42,.05);
}
.card:hover { border-color:var(--sea) }
.card .kicker { font-size:.72rem; font-weight:700; letter-spacing:.06em; text-transform:uppercase; color:var(--ember) }
.card h3 { font-family:Fraunces,serif; margin:.2rem 0 .35rem; font-size:1.18rem }
.card p { margin:0; font-size:.9rem }
.card .go { margin-top:.55rem; font-weight:600; color:var(--sea); font-size:.88rem }
.grid2 { display:grid; grid-template-columns:1fr 1fr; gap:.8rem }
.facts { display:grid; grid-template-columns:repeat(4,1fr); gap:.7rem; margin:0 0 1.15rem }
.fact { background:#fff; border:1px solid var(--line); border-radius:14px; padding:.85rem }
.fact b { display:block; font-family:Fraunces,serif; font-size:1.4rem; color:var(--deep) }
.fact span { color:var(--muted); font-size:.78rem }
.dim-table-wrap { overflow-x:auto; margin:0 0 1rem }
.dim-table {
  width:100%; border-collapse:collapse; font-size:.86rem; background:var(--card);
  border:1px solid var(--line); border-radius:12px; overflow:hidden;
}
.dim-table th, .dim-table td { padding:.55rem .7rem; text-align:left; border-bottom:1px solid var(--line) }
.dim-table th { background:rgba(15,118,110,.08); color:var(--deep); font-weight:600; font-size:.78rem }
.dim-table tr:last-child td { border-bottom:none }
.dim-table td.num { font-variant-numeric:tabular-nums; white-space:nowrap }
.dim-note { font-size:.82rem; color:var(--muted); max-width:64ch; margin:0 0 1rem }
.types { display:flex; flex-wrap:wrap; gap:.4rem; margin:.2rem 0 1.15rem }
.types span { background:#fff; border:1px solid var(--line); border-radius:999px; padding:.28rem .7rem; font-size:.78rem; color:var(--deep) }
.closing { margin-top:2.6rem; padding:1.4rem 1.35rem; border-radius:16px; background:linear-gradient(135deg,var(--deep),var(--ember)); color:#f5fffb }
.closing h2 { font-family:Fraunces,serif; margin:0 0 .55rem; font-size:1.4rem }
.closing p { margin:0; opacity:.95; max-width:62ch }
.footer { margin-top:1.6rem; color:var(--muted); font-size:.88rem }
.footer a { margin-right:.85rem }
.modal { display:none; position:fixed; inset:0; background:rgba(0,0,0,.88); z-index:30; align-items:center; justify-content:center; padding:1.2rem; cursor:zoom-out }
.modal.open { display:flex }
.modal img { max-width:95vw; max-height:90vh; border-radius:8px }
@media (max-width:760px) { .grid2, .facts { grid-template-columns:1fr 1fr } .hero { min-height:78vh } }
</style>
</head>
<body>
<header class="hero">
  <div>
    <p class="brand">Two ways to build a gap junction</p>
    <h1>A phylogenetic story of innexins and connexins</h1>
    <p class="lede">Trees, synteny, gene models, and a few places where the databases and GFF annotations get things wrong.</p>
    <div class="cta">
      <a href="#how">Start: how they work</a>
      <a class="ghost" href="#end">Skip to the contrast</a>
    </div>
  </div>
</header>

<nav class="path" aria-label="Story path">
  <div class="path-inner">
    <a href="#how">1. How they work</a>
    <a href="#inx">2. Innexins</a>
    <a href="#cnx">3. Connexins</a>
    <a href="#artifacts">4. Annotation issues</a>
    <a href="#inx-comp">5. Innexin vs innexin</a>
    <a href="#cnx-comp">6. Connexin vs connexin</a>
    <a href="#end">7. Innexin × connexin</a>
  </div>
</nav>

<main class="wrap">

  <section class="chapter" id="how">
    <span class="eyebrow">Chapter 1</span>
    <h2>Three families, one biological job</h2>
    <p>
      Neighbouring cells can exchange ions and small metabolites through gap junctions.
      In animals this is done by three protein families: innexins (mostly invertebrates),
      connexins (vertebrates), and pannexins (chordate/vertebrate innexin homologs that do not form
      classical intercellular junctions). This project compares them with alignments, trees, synteny,
      and gene models.
    </p>
    <div class="figure">
      <img src="figures/01_two_cells.png" alt="Two cells coupled by a gap junction" onclick="openModal(this.src)">
      <p class="cap">What the channels do in the cell. The rest of this site is about the genes behind them.</p>
    </div>
    <p>
      Innexins and connexins look similar in membrane topology (four transmembrane helices). Literature
      treats them as separate gap-junction channel families; this panel tests that with sequence markers
      and cross-family searches. A practical check is the conserved cysteine pattern in the extracellular
      loops: innexins and pannexins usually have 2+2 = 4, connexins 3+3 = 6. Pannexins also keep a
      conserved N-X-S/T motif on those loops, which marks the non-docking innexin-lineage branch.
    </p>
    <div class="dim-table-wrap">
      <table class="dim-table" aria-label="Sequence-family markers">
        <thead>
          <tr><th></th><th>Innexin</th><th>Connexin</th><th>Pannexin</th></tr>
        </thead>
        <tbody>
          <tr><td>Where they live</td><td class="num">Predominantly invertebrates</td><td class="num">Canonical vertebrate GJ family</td><td class="num">Chordates / vertebrates</td></tr>
          <tr><td>Cys motif (EL1+EL2)</td><td class="num">2 + 2 = 4</td><td class="num">3 + 3 = 6</td><td class="num">2 + 2 = 4</td></tr>
          <tr><td>Gene model (typical)</td><td class="num">More exons (median ~5)</td><td class="num">Compact ORF (~2 exons)</td><td class="num">Innexin-like</td></tr>
          <tr><td>Forms a gap junction?</td><td class="num">Yes</td><td class="num">Yes</td><td class="num">No (N-X-S/T marker)</td></tr>
        </tbody>
      </table>
    </div>
    <p class="dim-note">
      Structural details from cryo-EM are only background reading here. Main points on this path are
      the sequence patterns, the trees, and the annotation issues discussed in
      <a href="#artifacts">chapter 4</a>.
    </p>
    <div class="figure">
      <img src="figures/02_topology.png" alt="Four-transmembrane topology with 2 vs 3 cysteines per extracellular loop" onclick="openModal(this.src)">
      <p class="cap">Same rough fold idea, different Cys counts: 4 for innexin/pannexin, 6 for connexin.</p>
    </div>
    <div class="facts" aria-label="Conserved cysteine signature">
      <div class="fact"><b>2 + 2</b><span>inx / panx Cys</span></div>
      <div class="fact"><b>4</b><span>per innexin monomer</span></div>
      <div class="fact"><b>3 + 3</b><span>connexin Cys</span></div>
      <div class="fact"><b>6</b><span>per connexin monomer</span></div>
    </div>
    <div class="figure">
      <img src="figures/03_hexamers.png" alt="Octameric innexon versus hexameric connexon" onclick="openModal(this.src)">
      <p class="cap">Literature labels: innexon often drawn as an octamer, connexon as a hexamer.</p>
    </div>
    <div class="figure">
      <img src="figures/04_kingdoms.png" alt="Innexins in invertebrates, connexins in vertebrates" onclick="openModal(this.src)">
      <p class="cap">Kingdom occupancy in this panel: innexins in invertebrates, connexins in vertebrates — not one gene renamed twice.</p>
    </div>
    <p>
      <strong>Literature (Welzel &amp; Schuster 2022):</strong> lancelets apparently keep only one innexin with an
      N-X-S/T site and no connexins — an early chordate “innexin bottleneck” that may have left room for
      connexins to arise later, while the glycosylated innexin branch continued as pannexins.
      <strong>This panel:</strong> does not yet run a dedicated amphioxus copy-number / NGS audit; the sketch below
      is literature framing, not a result demonstrated here.
      <strong>Together:</strong> the bottleneck is a plausible bridge between innexin loss and connexin emergence,
      but the evidence on this site is occupancy + sequence separation, not a dated evolutionary reconstruction.
    </p>
    <div class="facts" aria-label="Chordate bottleneck">
      <div class="fact"><b>1</b><span>innexin in lancelets (lit.)</span></div>
      <div class="fact"><b>0</b><span>connexins in lancelets (lit.)</span></div>
      <div class="fact"><b>N-X-S/T</b><span>on that remnant</span></div>
      <div class="fact"><b>Cx</b><span>appear later (lit.)</span></div>
    </div>
    <div class="figure">
      <img src="figures/05_chordate_bottleneck.png" alt="Early chordate innexin bottleneck and connexin emergence" onclick="openModal(this.src)">
      <p class="cap">Sketch of the lancelet bottleneck idea from the literature.</p>
    </div>
    <a class="card" href="../../../pannexin/project/results/pannexin_path/index.html">
      <div class="kicker">Third family</div>
      <h3>Pannexin path</h3>
      <p>What pannexins are, PANX1–3, N-X-S/T, and links into the pannexin results.</p>
      <div class="go">Open →</div>
    </a>
    <a class="card" href="../../../pannexin/project/results/pannexin_clade_comparison/index.html" style="margin-top:.8rem">
      <div class="kicker">Clade × type</div>
      <h3>Pannexin clade comparison</h3>
      <p>PANX1–3 counts, copy number, lengths — same shape as connexin clade pages.</p>
      <div class="go">Open →</div>
    </a>
  </section>

  <section class="chapter" id="inx">
    <span class="eyebrow">Chapter 2</span>
    <h2>Innexins: names that do not equal orthology</h2>
    <p>
      Innexins predominate among invertebrate gap-junction channels in this panel (exceptions exist in
      chordates via pannexins). We searched 51 genomes and recovered
      142 present loci in 20 species (16 of those were new or previously blank). With the reference
      panel that is 164 loci in 28 species. Typical products are ~350–420 aa, four TM helices, and
      two conserved cysteines in each extracellular loop (4 per monomer). Coding-region gene models
      are intron-rich (panel median ~5 exons), so one locus can give several splice isoforms —
      familiar from Drosophila / mosquito shakB.
    </p>
    <div class="facts">
      <div class="fact"><b>51</b><span>genomes searched</span></div>
      <div class="fact"><b>20</b><span>species with present loci</span></div>
      <div class="fact"><b>142</b><span>present innexin loci</span></div>
      <div class="fact"><b>4</b><span>tree subfamilies (SF1–SF4)</span></div>
    </div>
    <p>
      Named types in the panel are the Drosophila set Inx1/ogre, Inx2, Inx3, Inx4 (zpg), Inx5, Inx6, Inx7, shakB
      and the nematode set unc-7, unc-9, inx-2 (plus other nematode genes).
      Those names are labels from classical screens, not orthology claims: recovered non-insect loci
      most often best-hit unc-9 (39) or unc-7 (34). Orthology here is read from the tree and from synteny.
    </p>
    <div class="types" aria-label="Innexin types">
      <span>Inx1 / ogre</span><span>Inx2</span><span>Inx3</span><span>Inx4 / zpg</span><span>Inx5</span><span>Inx6</span><span>Inx7</span><span>shakB</span><span>unc-7</span><span>unc-9</span><span>inx-2 (Cele)</span>
    </div>
    <p>
      The tree does not replay “Inx1–7 + shakB” in every animal. Four subfamilies cut across gene names:
      SF1 shakB, SF2 ogre/Inx1, SF3 Inx3 + Inx7 + all nematode innexins, SF4 Inx2/zpg/Inx4–6.
      Outside insects, recovered loci look like the SF3 / nematode-like radiation.
      Drosophila even places Inx7 next to ogre/Inx2 on the X — genomic neighbours, phylogenetic strangers.
      Insect Inx2 (SF4) and nematode inx-2 (SF3) share a numerical suffix from screens but are not orthologs;
      the useful contrast is tree-guided SF3↔nematode synteny
      (<a href="../innexin_subfamilies/SF3_Inx3_Inx7_vs_nematode_innexins.html">SF3 panel</a>), not name matching.
    </p>
    <a class="card" href="../innexin_insights/index.html">
      <div class="kicker">Overview</div>
      <h3>Innexin insights</h3>
      <p>Short read of the innexin tree: clade geography, identity, Drosophila architecture, embeddings.</p>
      <div class="go">Open →</div>
    </a>
    <div class="cards grid2">
      <a class="card" href="../new_species_gallery/index.html">
        <div class="kicker">Recoveries</div>
        <h3>New species gallery</h3>
        <p>16 new/rescued species, 142 present loci.</p>
        <div class="go">Open →</div>
      </a>
      <a class="card" href="../innexin_subfamilies/index.html">
        <div class="kicker">Tree labels</div>
        <h3>Four subfamilies</h3>
        <p>SF1–SF4, and why nematode inx-2 is not insect Inx2.</p>
        <div class="go">Open →</div>
      </a>
    </div>
  </section>

  <section class="chapter" id="cnx">
    <span class="eyebrow">Chapter 3</span>
    <h2>Connexins: named classes in stable neighbourhoods</h2>
    <p>
      Connexins are the vertebrate gap-junction proteins. In the literature they appear after the early
      chordate innexin bottleneck (lancelets still lack them). This panel has 103 reference proteins
      plus 24 discovery loci across 38 species. Typical products are compact (~250–380 aa), four TM helices,
      and three conserved cysteines in each extracellular loop (6 per monomer).
      Unlike innexins, the coding region usually has no introns (ORF in one exon); models often still
      report ~2 exons because of a separate 5′ UTR exon (panel median 2).
    </p>
    <div class="facts">
      <div class="fact"><b>38</b><span>vertebrate species</span></div>
      <div class="fact"><b>103</b><span>reference connexins</span></div>
      <div class="fact"><b>24</b><span>discovery loci</span></div>
      <div class="fact"><b>5</b><span>classes (α β γ δ ε)</span></div>
    </div>
    <p>
      Named types follow GJA/GJB/GJC/GJD/GJE:
      GJA1, GJA3, GJA4, GJA5, GJA8 (α); GJB1, GJB2, GJB3, GJB6 (β); plus GJC, GJD, GJE.
      GJA1 is the most widely shared α gene; mammals are the densest part of the panel.
      Untyped discovery sequences stay unassigned — they are not a hidden innexin-like radiation.
    </p>
    <div class="types" aria-label="Connexin types">
      <span>GJA1</span><span>GJA3</span><span>GJA4</span><span>GJA5</span><span>GJA8</span><span>GJB1</span><span>GJB2</span><span>GJB3</span><span>GJB6</span><span>GJC</span><span>GJD</span><span>GJE</span>
    </div>
    <p>
      α (GJA), β (GJB), γ (GJC) and δ (GJD) co-occur across mammals, birds, amphibians and fish.
      GJA4 sits in a conserved β-cluster (GJB3–GJB5, DLGAP3, SMIM12) — a tandem-duplication array from
      local unequal crossing-over — while whole-genome duplications also put other paralogs on separate
      chromosomes (notably GJA1). Between <em>Xenopus laevis</em> (allotetraploid, L) and <em>X. tropicalis</em>
      (diploid) there is a chromosomal inversion around gja4: in laevis the flank
      <em>kdm1a.L → tmem30b.L → zmpste24.L → smim12.L</em> is left of <em>gja4.L</em>; in tropicalis that block
      is inverted and sits to the right of <em>gja4</em>. On the
      <a href="../synteny_gff/cross_species/connexin_GJA1_cross_species.html">GJA1 microsynteny</a> plot,
      empty flanks around <em>X. laevis</em> <em>gja1.L</em> look like an isolated gene but are a short scaffold
      or assembly gap — compare with the rich Danio neighbourhood.
      Within-class identity is ~57% vs ~40% between; α↔β ~46%.
    </p>
    <a class="card" href="../connexin_insights/index.html">
      <div class="kicker">Overview</div>
      <h3>Connexin insights</h3>
      <p>α/β co-occurrence, GJA4 neighbourhoods, embeddings, compact exons.</p>
      <div class="go">Open →</div>
    </a>
    <div class="cards grid2">
      <a class="card" href="../synteny_gff/cross_species/connexin_GJA4_cross_species.html">
        <div class="kicker">Locus</div>
        <h3>GJA4 microsynteny</h3>
        <p>Tandem β-cluster and the Xenopus inversion.</p>
        <div class="go">Open →</div>
      </a>
      <a class="card" href="../synteny_gff/cross_species/connexin_GJA1_cross_species.html">
        <div class="kicker">Locus</div>
        <h3>GJA1 microsynteny</h3>
        <p>WGD-dispersed α address; empty Xenopus flanks are a scaffold/gap issue.</p>
        <div class="go">Open →</div>
      </a>
    </div>
  </section>

  <section class="chapter" id="artifacts">
    <span class="eyebrow">Chapter 4</span>
    <h2>Where the panels look wrong</h2>
    <p>
      A few numbers in the plots look like biology until you check the gene models and the sampling.
      Three of them are annotation or panel artifacts. Two synteny examples below are real rearrangements.
    </p>

    <h3 style="font-family:Fraunces,serif;font-size:1.15rem;margin:1.2rem 0 .45rem;color:var(--deep)">1. δ GJD “~22 exons”</h3>
    <p>
      The connexin exon-by-subfamily boxplot can spike toward ~22 exons in δ (GJD). That is not a real
      22-exon connexin. Connexin ORFs are usually intron-free (~2 exons once a 5′ UTR exon is counted).
      Automated GFF often fuses that compact ORF with intron-rich neighbours in the GJA4 neighbourhood
      (DLGAP3, SMIM12). Extreme exon counts need a manual look at the locus.
      Plot: <a href="../connexin_clade_comparison/index.html">connexin clade comparison</a>.
    </p>

    <h3 style="font-family:Fraunces,serif;font-size:1.15rem;margin:1.2rem 0 .45rem;color:var(--deep)">2. “1-exon” innexin median (Ctenophora / non-insect Arthropoda)</h3>
    <p>
      Outside insects the exon-by-clade plot can show a median of exactly 1.0. Real innexins (insects,
      nematodes) have coding-region introns and splice isoforms (panel median ~5 exons). The flat 1-exon
      median mostly comes from TSA or ab initio models that report one continuous ORF when intron
      boundaries were never resolved.
      Plot: <a href="../innexin_clade_comparison/index.html">innexin clade comparison</a>.
    </p>

    <h3 style="font-family:Fraunces,serif;font-size:1.15rem;margin:1.2rem 0 .45rem;color:var(--deep)">3. Connexin copy number ≈ 1</h3>
    <p>
      An early histogram gave a median of ~1 connexin locus per species. That was the UniProt / discovery
      slice we used (23/38 species had only one locus in that panel), not the full vertebrate repertoire.
      Jawed vertebrates typically have ~20–22 connexin genes (human 21; our panel human 22, mouse 20);
      teleosts can reach ~46 after extra WGDs. The corrected figure shows well-covered reference species
      against that literature band.
      Plot: <a href="../inx_vs_cnx/index.html#copy">inx × cnx · copy number</a>.
    </p>

    <h3 style="font-family:Fraunces,serif;font-size:1.15rem;margin:1.2rem 0 .45rem;color:var(--deep)">4. Other plot traps</h3>
    <ul style="color:var(--muted);max-width:68ch;margin:0 0 1rem;padding-left:1.15rem">
      <li style="margin-bottom:.45rem">
        Empty flanks around <em>X. laevis</em> <em>gja1.L</em> vs a rich Danio neighbourhood: short scaffold
        or assembly gap, not a biological gene island.
        <a href="../synteny_gff/cross_species/connexin_GJA1_cross_species.html">GJA1 microsynteny</a>.
      </li>
      <li style="margin-bottom:.45rem">
        Insect Inx2 (SF4) ≠ nematode inx-2 (SF3): same numerical suffix from screens, different clades.
        Use SF3↔nematode synteny, not the names.
        <a href="../innexin_subfamilies/SF3_Inx3_Inx7_vs_nematode_innexins.html">SF3 panel</a>.
      </li>
    </ul>

    <h3 style="font-family:Fraunces,serif;font-size:1.15rem;margin:1.2rem 0 .45rem;color:var(--deep)">5. Synteny that is biology</h3>
    <ul style="color:var(--muted);max-width:68ch;margin:0 0 1rem;padding-left:1.15rem">
      <li style="margin-bottom:.45rem">
        Xenopus GJA4 inversion: allotetraploid <em>X. laevis</em> has
        <em>kdm1a.L → tmem30b.L → zmpste24.L → smim12.L</em> left of <em>gja4.L</em>; diploid
        <em>X. tropicalis</em> inverts that block to the right of <em>gja4</em>. Same neighbourhood,
        different orientation. The panel is also a tandem array (cx38 / gjb3 / gjb4 / gja4), unlike
        WGD-dispersed GJA1.
        <a href="../synteny_gff/cross_species/connexin_GJA4_cross_species.html">GJA4 microsynteny</a>.
      </li>
      <li style="margin-bottom:.45rem">
        Drosophila Inx7 sits next to ogre/Inx2 on the X, but the tree puts Inx7 with nematode innexins
        in SF3. Cross-phylum synteny should compare Dmel Inx3/Inx7 (SF3) to nematode innexins as a group,
        not Inx2 vs inx-2 by name.
      </li>
    </ul>

    <div class="facts" aria-label="Annotation caveats">
      <div class="fact"><b>22 exons</b><span>GFF fusion · δ GJD</span></div>
      <div class="fact"><b>1 exon</b><span>TSA / ab initio</span></div>
      <div class="fact"><b>median 1</b><span>copy-number slice</span></div>
      <div class="fact"><b>empty GJA1</b><span>scaffold / gap</span></div>
    </div>
  </section>

  <section class="chapter" id="inx-comp">
    <span class="eyebrow">Chapter 5</span>
    <h2>Innexin subfamilies in numbers</h2>
    <p>
      Within-subfamily identity is ~46% versus ~28% between subfamilies. SF3↔SF4 sits at ~28% —
      insect expansion versus nematode-like radiation. 3-mer embeddings agree (same group ~0.28 vs
      different ~0.05); amino-acid composition does not. Outside insects, clade × best-hit maps pile
      into unc-7 / unc-9 / SF3, not a full Drosophila toolkit.
    </p>
    <a class="card" href="../innexin_insights/index.html">
      <div class="kicker">Overview</div>
      <h3>Innexin insights</h3>
      <p>Clade geography, identity, embeddings, length and exons.</p>
      <div class="go">Open →</div>
    </a>
    <div class="cards grid2">
      <a class="card" href="../innexin_clade_comparison/index.html">
        <div class="kicker">Tables</div>
        <h3>Innexin clade comparison</h3>
        <p>Every locus by clade and best-hit label.</p>
        <div class="go">Open →</div>
      </a>
      <a class="card" href="../innexin_similarity/index.html">
        <div class="kicker">Similarity</div>
        <h3>Innexin MMseqs + embeddings</h3>
        <p>All-vs-all identity and 3-mer PCA (46% vs 28%).</p>
        <div class="go">Open →</div>
      </a>
    </div>
  </section>

  <section class="chapter" id="cnx-comp">
    <span class="eyebrow">Chapter 6</span>
    <h2>Connexin classes in numbers</h2>
    <p>
      Within-class identity ~57% versus ~40% between classes; α↔β ~46%. Embeddings split α from β
      more cleanly than composition. Named GJA/GJB labels transfer across vertebrate clades —
      unlike insect Inx1–7 names that do not transfer to nematodes.
    </p>
    <a class="card" href="../connexin_insights/index.html">
      <div class="kicker">Overview</div>
      <h3>Connexin insights</h3>
      <p>α/β co-occurrence, identity heatmaps, embeddings, gene models.</p>
      <div class="go">Open →</div>
    </a>
    <div class="cards grid2">
      <a class="card" href="../connexin_clade_comparison/index.html">
        <div class="kicker">Tables</div>
        <h3>Connexin clade comparison</h3>
        <p>Which GJA/GJB/GJC labels appear in which vertebrate clade.</p>
        <div class="go">Open →</div>
      </a>
      <a class="card" href="../connexin_similarity/index.html">
        <div class="kicker">Similarity</div>
        <h3>Connexin MMseqs + embeddings</h3>
        <p>57% vs 40%, and why α↔β is still one family.</p>
        <div class="go">Open →</div>
      </a>
    </div>
  </section>

  <section class="chapter" id="end">
    <span class="eyebrow">Chapter 7</span>
    <h2>Innexin × connexin</h2>
    <p>
      Same axes for both families. Different kingdoms. In one 3-mer space they form two clouds
      (innexin↔connexin cosine −0.40). MMseqs finds no homologous alignments. Lengths differ
      (median 367 vs 322 aa). Gene models differ: innexins carry coding-region introns and splice
      isoforms; connexin ORFs are usually intron-free. Genomes differ too: tandem clusters that cut
      the innexin tree, versus conserved connexin neighbourhoods.
    </p>
    <div class="figure">
      <img src="../inx_vs_cnx/figures/02_joint_kmer_pca.png" alt="Joint 3-mer PCA of innexins and connexins" onclick="openModal(this.src)">
      <p class="cap">Both families in one sequence space (from the direct-contrast page).</p>
    </div>
    <a class="card" href="../inx_vs_cnx/index.html">
      <div class="kicker">Main contrast</div>
      <h3>Direct innexin × connexin</h3>
      <p>Kingdom occupancy, embeddings, homology, length, exons, copy number, neighbourhoods.</p>
      <div class="go">Open →</div>
    </a>
    <a class="card" href="../family_comparison/index.html" style="margin-top:.8rem">
      <div class="kicker">Appendix</div>
      <h3>Same metrics, side by side</h3>
      <p>Within/between identity for each family on one page — not the direct contrast.</p>
      <div class="go">Open →</div>
    </a>
  </section>

  <section class="closing">
    <h2>Summary</h2>
    <p>
      <strong>Literature:</strong> innexin- and connexin-based gap junctions are often treated as separate
      evolutionary solutions, with an early chordate innexin bottleneck (Welzel &amp; Schuster 2022) preceding
      connexin diversification and a retained glycosylated innexin branch (pannexins).
      <strong>This panel:</strong> innexins diversify mainly in invertebrates; connexins appear as named vertebrate
      classes in stable neighbourhoods; cross-family MMseqs finds 0 hits at e≤1e-3 and ≥80 aa and the two
      families occupy separate clouds in sequence space.
      <strong>Together:</strong> parallel gap-junction toolkits — not one renamed homologous family — but the dated
      bottleneck → pannexin → connexin sequence is literature context, not something re-derived from these trees alone.
    </p>
  </section>

  <footer class="footer">
    <p>Reading path over the result pages.</p>
    <p>
      <a href="tools.html">Tools used</a>
      <a href="../phylogenetic_story/index.html">Phylogeny story</a>
      <a href="../showcase/index.html">Full showcase</a>
      <a href="../../../pannexin/project/results/pannexin_path/index.html">Pannexin path</a>
      <a href="../../../pannexin/project/results/panx_vs_inx/index.html">Pannexin × innexin</a>
      <a href="../sources/index.html">Sources</a>
    </p>
  </footer>
</main>

<div id="modal" class="modal" onclick="this.classList.remove('open')"><img id="modal-img" alt=""></div>
<script>
function openModal(src) {
  document.getElementById('modal-img').src = src;
  document.getElementById('modal').classList.add('open');
}
document.addEventListener('keydown', e => {
  if (e.key === 'Escape') document.getElementById('modal').classList.remove('open');
});
</script>
</body>
</html>
'''
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "index.html").write_text(html, encoding="utf-8")


def write_tools_html() -> None:
    """Methods note: tools and parameters used in this project."""
    html = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Tools &amp; parameters — gap-junction analyses</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,560;9..144,700&family=Sora:wght@400;500;600&display=swap" rel="stylesheet">
<style>
:root {
  --ink:#12202a; --muted:#4d6570; --deep:#0a3a42; --sea:#0f766e; --ember:#c2410c;
  --foam:#f4f7f5; --card:#fff; --line:rgba(18,32,42,.1);
}
* { box-sizing:border-box }
body {
  margin:0; color:var(--ink); font-family:Sora,sans-serif; line-height:1.6;
  background:linear-gradient(180deg,#e8f1ee 0%, var(--foam) 40%, #eef3f1 100%);
}
a { color:var(--sea) }
.hero {
  padding:clamp(1.4rem,4vw,2.4rem); color:#f7fffb;
  background:linear-gradient(118deg, rgba(10,58,66,.94) 0%, rgba(15,118,110,.72) 55%, rgba(194,65,12,.55) 100%);
}
.hero p { font-family:Fraunces,serif; font-size:clamp(1.8rem,4vw,2.6rem); margin:0 0 .4rem; letter-spacing:-.03em }
.hero span { opacity:.9; font-size:.95rem }
.wrap { width:min(780px, calc(100% - 2rem)); margin:0 auto; padding:1.6rem 0 3.5rem }
h2 { font-family:Fraunces,serif; font-size:1.15rem; margin:1.7rem 0 .55rem }
ul { margin:0; padding:0; list-style:none }
li {
  background:var(--card); border:1px solid var(--line); border-radius:12px;
  padding:.75rem .95rem; margin:0 0 .45rem; font-size:.9rem;
}
li b { color:var(--deep) }
li span { display:block; color:var(--muted); font-size:.82rem; margin-top:.2rem }
code {
  font-family:ui-monospace,SFMono-Regular,Menlo,monospace; font-size:.8em;
  background:rgba(15,118,110,.08); padding:.08em .28em; border-radius:4px;
}
.params {
  margin-top:.45rem; padding:.45rem .55rem; border-left:3px solid rgba(15,118,110,.35);
  background:rgba(15,118,110,.04); border-radius:0 8px 8px 0; font-size:.8rem; color:var(--ink);
  line-height:1.55;
}
.back { display:inline-block; margin:.4rem 0 0; color:#f7fffb; font-weight:600; text-decoration:none }
.note { color:var(--muted); font-size:.85rem; margin-top:1.4rem }
</style>
</head>
<body>
<header class="hero">
  <p>Tools &amp; parameters</p>
  <span>Software and settings used for the diagrams and analyses on this path.</span>
  <div><a class="back" href="index.html">← back to the story</a></div>
</header>
<main class="wrap">

  <h2>Finding genes in genomes</h2>
  <ul>
    <li>
      <b>miniprot</b>
      <span>Coarse location of innexin/connexin-like loci on genome assemblies (discovery + curator search).</span>
      <div class="params">
        <strong>Discovery</strong> (<code>discover_innexins.py</code>): <code>-j 2</code> (insect/vertebrate splice model);
        score floor ≥ 30; HSP filters ≥ 60 bp, ≥ 20 aa query, ≥ 18% identity (short hits ≥ 35 aa and ≥ 22%);
        locus clustering gap ≤ 8 kb with ±2 kb flanks; large genomes (&gt;4 GB) run per-contig with
        contig ≥ 1 Mb and up to 3 parallel jobs.<br>
        <strong>Curator probe</strong>: <code>--outs=0.5</code>, 4 threads; keep models with identity ≥ 0.20
        and length ≥ 280 aa / ≥ 4 Cys (present), ≥ 180 aa (fragmentary), or ≥ 100 aa (weak).
      </div>
    </li>
    <li>
      <b>samtools</b>
      <span>FAI indexing and region extraction for translating miniprot models.</span>
      <div class="params"><code>samtools faidx</code> on genome FASTAs; region pull for CDS spans before Biopython translation.</div>
    </li>
    <li>
      <b>MMseqs2</b> — candidate typing
      <span>Protein search of candidate ORFs against trusted reference innexins/connexins.</span>
      <div class="params">
        <code>easy-search</code>, <code>--search-type 1</code> (protein), <code>-e 1e-5</code>,
        <code>-s 5.7</code>, format <code>query,target,pident,qcov,bits</code>.
        Keep ranks with reference identity ≥ 15%; protein length window ~150–750 aa (innexin) / similar for connexin.
      </div>
    </li>
  </ul>

  <h2>Alignment, trees, similarity</h2>
  <ul>
    <li>
      <b>MAFFT</b>
      <span>Multiple sequence alignment of phylogeny input FASTAs.</span>
      <div class="params">Default mode <code>mafft --auto</code> (optional <code>linsi</code> / <code>ginsi</code> / <code>fftns</code>).</div>
    </li>
    <li>
      <b>IQ-TREE 2</b>
      <span>Phylogenetic trees for the innexin and connexin stories.</span>
      <div class="params">
        <code>-m MFP</code> (ModelFinder Plus); <code>-bb 1000</code> ultrafast bootstrap;
        <code>-nt 2</code> threads (default in runner).
      </div>
    </li>
    <li>
      <b>MMseqs2</b> — all-vs-all &amp; clustering
      <span>Within-family identity heatmaps and threshold sweeps (innexin / connexin similarity pages).</span>
      <div class="params">
        <strong>Search:</strong> <code>mmseqs search</code> all-vs-all, <code>-s 7.5</code>, <code>-e 1e-3</code>,
        <code>--max-seqs 1000</code> (innexin) / <code>800</code> (connexin).<br>
        <strong>Cluster sweep:</strong> <code>--min-seq-id</code>
        innexin <code>0.20–0.70</code> (step 0.10); connexin <code>0.30–0.80</code>;
        coverage <code>-c 0.8</code>, <code>--cov-mode 0</code>.
      </div>
    </li>
    <li>
      <b>MMseqs2</b> — cross-family &amp; discovery typing
      <span>Innexin↔connexin homology check; discovery loci typed against reference connexins.</span>
      <div class="params">
        Cross-family: <code>easy-search</code> innexin query → connexin target, <code>-e 10</code>,
        <code>-s 7.5</code>, <code>--max-seqs 5</code>. Homologous call =
        <strong>e-value ≤ 1e-3 and alignment length ≥ 80 aa</strong> (raw short fragments are not counted as homology).<br>
        Discovery typing: <code>easy-search</code>, <code>-e 1e-3</code>, <code>--max-seqs 5</code>.
      </div>
    </li>
    <li>
      <b>NumPy embeddings</b>
      <span>3-mer composition + SVD and amino-acid composition cosines (family contrast pages).</span>
      <div class="params">
        k-mers: <code>k=3</code>, vocabulary capped at top 800 kmers, truncated SVD to
        <code>n_comp=20</code> (similarity pages) or <code>8</code> (direct innexin×connexin PCA);
        vectors L2-normalised; mean pairwise cosine sampled (seed 7, ≤ 8000 pairs).
      </div>
    </li>
  </ul>

  <h2>Synteny and gene architecture</h2>
  <ul>
    <li>
      <b>SynVoy</b>
      <span>Synteny-guided neighbourhood comparisons around selected innexin and connexin loci.</span>
      <div class="params">
        Easy mode: <code>--n_flanking_genes 5</code>, <code>--adaptive_max_regions 3</code>,
        <code>--enable_smith_waterman false</code>, <code>--exon_level_search false</code>,
        <code>--auto_params false</code>, <code>--multi_profile false</code>,
        <code>--max_genomes 3</code> (local queue). Reference insect runs also use
        <code>--mmseqs_split_memory_limit 1G</code>.
      </div>
    </li>
    <li>
      <b>NCBI GFF microsynteny</b>
      <span>Exon counts, gene span, and cross-species neighbour sharing (GJA1 / GJA4 / innexin panels).</span>
      <div class="params">
        Neighbour plots: ±180 kb around the gene of interest, ≤ 22 genes kept in the window.
        Neighbour conservation tiers from shared named genes across species (high = ≥3 shared in all species).
      </div>
    </li>
    <li>
      <b>UniProt / Pfam</b>
      <span>Reference protein sequences and domain labels (e.g. Connexin PF00029).</span>
    </li>
  </ul>

  <h2>Family sequence markers</h2>
  <ul>
    <li>
      <b>Extracellular cysteine motif</b>
      <span>Conserved Cys pattern in the two extracellular loops (alignment filter).</span>
      <div class="params">
        Innexin / pannexin: 2 Cys per loop → 4 per monomer.<br>
        Connexin: 3 Cys per loop → 6 per monomer.<br>
        Filters that require ≥4 Cys are a product-cleanliness check. Whole-protein Cys totals are often higher.
      </div>
    </li>
    <li>
      <b>Hemichannel stoichiometry (literature)</b>
      <span>Not measured in this panel — labels used when reading papers.</span>
      <div class="params">
        Innexon = octamer (8). Connexon = hexamer (6). Pannexon = heptamer (7; non-docking).
      </div>
    </li>
    <li>
      <b>Coding-region introns / splicing</b>
      <span>From reference GFF gene models (<code>exon_structures/</code>).</span>
      <div class="params">
        Connexin: ORF usually uninterrupted; models often ~2 exons (5′ UTR + coding exon). Panel median = 2.<br>
        Innexin: coding-region introns common (panel median ~5 exons); one locus can yield several isoforms (e.g. <em>shakB</em>).
        GFF exon counts include UTR exons where annotated.
      </div>
    </li>
    <li>
      <b>Annotation issues in the plots</b>
      <span>See also <a href="index.html#artifacts">chapter 4</a> on the story path.</span>
      <div class="params">
        Single-exon innexin medians (Ctenophora / non-insect Arthropoda): often TSA / ab initio ORFs without resolved introns.<br>
        δ GJD “~22 exons”: GFF fusion of a compact connexin ORF with intron-rich flanks (DLGAP3 / SMIM12 near GJA4).<br>
        Empty GJA1 flanks in <em>X. laevis</em>: short scaffold / assembly gap, not a gene island.<br>
        Insect Inx2 (SF4) ≠ nematode inx-2 (SF3): same name, different clades — use tree-guided SF3 contrasts.
      </div>
    </li>
    <li>
      <b>Early chordate bottleneck (literature)</b>
      <span>Welzel &amp; Schuster 2022 (<em>eLife</em>) — gene-content framing, not a panel parameter.</span>
      <div class="params">
        Lancelets lack connexins and keep one innexin with an extracellular N-X-S/T motif (pannexin-lineage marker).
        That bottleneck is the usual explanation for later connexin diversification.
        See <a href="../../../pannexin/project/results/pannexin_path/index.html">pannexin path</a>.
      </div>
    </li>
  </ul>

  <h2>Figures and pages</h2>
  <ul>
    <li><b>Python · Matplotlib</b><span>Custom diagrams on this path. Agg backend for batch export.</span></li>
    <li><b>pandas · Biopython</b><span>Tables, FASTA I/O, and sequence translation from miniprot CDS models.</span></li>
  </ul>

  <p class="note">Data sources: NCBI genome assemblies and UniProt reference FASTAs. Reported % identity values are MMseqs2 pairwise identities, not BLAST.</p>
  <p class="note"><a href="index.html">← story path</a></p>
</main>
</body>
</html>
'''
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "tools.html").write_text(html, encoding="utf-8")


def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    plot_two_cells(FIGS / "01_two_cells.png")
    plot_topology(FIGS / "02_topology.png")
    plot_hexamers(FIGS / "03_hexamers.png")
    plot_kingdoms(FIGS / "04_kingdoms.png")
    plot_chordate_bottleneck(FIGS / "05_chordate_bottleneck.png")
    plot_more_than_pipes(FIGS / "06_more_than_pipes.png")
    plot_biophysical_contrasts(FIGS / "07_biophysical_contrasts.png")
    write_html()
    write_tools_html()
    print(f"Done → {OUT / 'index.html'}")
    print(f"Tools → {OUT / 'tools.html'}")


if __name__ == "__main__":
    main()
