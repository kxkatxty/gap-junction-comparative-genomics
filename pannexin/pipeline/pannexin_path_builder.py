#!/usr/bin/env python3
"""Guided reading path for the pannexin project (separate from gap-junction thesis)."""

from __future__ import annotations

from pathlib import Path
import shutil

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import Circle, FancyBboxPatch, Rectangle

from pipeline.common import RESULTS_DIR

OUT = RESULTS_DIR / "pannexin_path"
FIGS = OUT / "figures"


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_channel_vs_junction(path: Path) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(9.2, 3.6))
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 6)
        ax.axis("off")

    # Left: pannexin hemichannel
    ax = axes[0]
    ax.add_patch(Rectangle((0, 2.2), 10, 1.6, color="#dbe4ea"))
    ax.text(5, 5.3, "Pannexin channel", ha="center", fontsize=12, fontweight="bold", color="#1e3a4c")
    ax.add_patch(Circle((5, 3.0), 0.85, facecolor="#c2410c", edgecolor="#7c2d12", lw=1.4, alpha=0.9))
    ax.add_patch(Circle((5, 3.0), 0.28, facecolor="#f8fafc", edgecolor="#7c2d12", lw=1))
    # N-glycan stubs on the extracellular face
    for dx in (-0.55, 0.55):
        ax.add_patch(Circle((5 + dx, 3.95), 0.18, facecolor="#fbbf24", edgecolor="#b45309", lw=1, zorder=5))
        ax.add_patch(Circle((5 + dx, 4.25), 0.12, facecolor="#fde68a", edgecolor="#b45309", lw=0.8, zorder=5))
    ax.text(5, 4.85, "N-glycans on EL", ha="center", fontsize=8, color="#b45309")
    ax.annotate("ATP / ions", xy=(5, 3.55), xytext=(7.6, 4.5),
                arrowprops=dict(arrowstyle="->", color="#475569"), fontsize=9, color="#475569")
    ax.text(5, 0.7, "one membrane → extracellular space", ha="center", fontsize=9, color="#64748b")

    # Right: classical gap junction (for contrast)
    ax = axes[1]
    ax.add_patch(Rectangle((0, 3.1), 10, 1.2, color="#dbe4ea"))
    ax.add_patch(Rectangle((0, 1.7), 10, 1.2, color="#e8efe9"))
    ax.text(5, 5.3, "Classical gap junction", ha="center", fontsize=12, fontweight="bold", color="#1e3a4c")
    ax.add_patch(Circle((5, 3.7), 0.7, facecolor="#0f766e", edgecolor="#115e59", lw=1.4, alpha=0.9))
    ax.add_patch(Circle((5, 2.3), 0.7, facecolor="#0f766e", edgecolor="#115e59", lw=1.4, alpha=0.9))
    ax.add_patch(Circle((5, 3.0), 0.22, facecolor="#f8fafc", edgecolor="#115e59", lw=1))
    ax.text(5, 0.7, "two cells coupled (innexin / connexin job)", ha="center", fontsize=9, color="#64748b")
    _save(fig, path)


def plot_glycosylation_blockade(path: Path) -> None:
    """N-linked glycosylation on extracellular loops blocks intercellular docking."""
    fig, axes = plt.subplots(1, 2, figsize=(9.6, 4.2))
    for ax in axes:
        ax.set_xlim(0, 10)
        ax.set_ylim(0, 7.2)
        ax.axis("off")

    def membrane_pair(ax, y_top: float = 4.35, y_bot: float = 1.85) -> None:
        ax.add_patch(Rectangle((0, y_top), 10, 1.2, color="#dbe4ea", zorder=0))
        ax.add_patch(Rectangle((0, y_bot), 10, 1.2, color="#e8efe9", zorder=0))
        ax.text(0.4, y_top + 0.6, "cell A", fontsize=8.5, color="#475569", va="center", fontweight="bold")
        ax.text(0.4, y_bot + 0.6, "cell B", fontsize=8.5, color="#475569", va="center", fontweight="bold")

    def glycan_cluster(ax, x: float, y: float) -> None:
        """Yellow N-glycan stub on an extracellular loop."""
        ax.add_patch(Circle((x, y), 0.2, facecolor="#fbbf24", edgecolor="#b45309", lw=1.1, zorder=6))
        ax.add_patch(Circle((x - 0.15, y + 0.22), 0.12, facecolor="#fde68a", edgecolor="#b45309", lw=0.8, zorder=6))
        ax.add_patch(Circle((x + 0.15, y + 0.22), 0.12, facecolor="#fde68a", edgecolor="#b45309", lw=0.8, zorder=6))

    label_kw = dict(
        ha="center",
        va="center",
        fontsize=10,
        fontweight="bold",
        zorder=12,
        clip_on=False,
    )

    # Left: two opposed pannexons — NGS on both extracellular docking faces
    ax = axes[0]
    membrane_pair(ax)
    ax.text(5, 6.55, "Pannexins — docking blocked", ha="center", fontsize=11.5, fontweight="bold", color="#c2410c")
    ax.add_patch(Circle((5, 4.9), 0.7, facecolor="#c2410c", edgecolor="#7c2d12", lw=1.5, alpha=0.92, zorder=3))
    ax.add_patch(Circle((5, 2.5), 0.7, facecolor="#c2410c", edgecolor="#7c2d12", lw=1.5, alpha=0.92, zorder=3))
    ax.add_patch(Circle((5, 4.9), 0.22, facecolor="#f8fafc", edgecolor="#7c2d12", lw=1, zorder=4))
    ax.add_patch(Circle((5, 2.5), 0.22, facecolor="#f8fafc", edgecolor="#7c2d12", lw=1, zorder=4))
    # Glycans face into the intercellular gap (EL of each hemichannel) — both sides
    for dx in (-0.42, 0.42):
        glycan_cluster(ax, 5 + dx, 4.2)   # cell A EL → gap
        glycan_cluster(ax, 5 + dx, 3.2)   # cell B EL → gap
    ax.text(5, 3.7, "N-glycans", ha="center", fontsize=8, color="#b45309", zorder=7)
    ax.text(
        7.55,
        3.55,
        "✕  no dock",
        color="#b45309",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#fff7ed", edgecolor="#fdba74", linewidth=1.2),
        **label_kw,
    )
    ax.text(
        5,
        0.55,
        "NGS on extracellular loops of both hemichannels\nsterically blocks head-to-head docking",
        ha="center",
        fontsize=8.5,
        color="#475569",
    )

    # Right: docked innexin/connexin GJ — hemichannels meet directly (no connector bar)
    ax = axes[1]
    membrane_pair(ax)
    ax.text(5, 6.55, "Innexin / connexin — docked GJ", ha="center", fontsize=11.5, fontweight="bold", color="#0f766e")
    ax.add_patch(Circle((5, 4.55), 0.72, facecolor="#0f766e", edgecolor="#115e59", lw=1.5, alpha=0.92, zorder=3))
    ax.add_patch(Circle((5, 2.85), 0.72, facecolor="#0f766e", edgecolor="#115e59", lw=1.5, alpha=0.92, zorder=3))
    ax.add_patch(Circle((5, 3.7), 0.22, facecolor="#f8fafc", edgecolor="#115e59", lw=1, zorder=4))
    ax.text(
        7.55,
        3.7,
        "✓  docked",
        color="#115e59",
        bbox=dict(boxstyle="round,pad=0.35", facecolor="#ecfdf5", edgecolor="#6ee7b7", linewidth=1.2),
        **label_kw,
    )
    ax.text(
        5,
        0.55,
        "No blocking NGS on the docking face —\nhemichannels couple across the gap",
        ha="center",
        fontsize=8.5,
        color="#475569",
    )
    _save(fig, path)


def plot_three_paralogs(path: Path) -> None:
    fig, ax = plt.subplots(figsize=(8.4, 3.4))
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.text(6, 4.4, "Vertebrate pannexin paralogs", ha="center", fontsize=13, fontweight="bold", color="#1e3a4c")
    cards = [
        (1.2, "PANX1", "#C2410C", "broad expression\nATP release / inflammation"),
        (4.7, "PANX2", "#1D4ED8", "more divergent\nstrong in nervous system"),
        (8.2, "PANX3", "#15803D", "closer to PANX1\nbone / skin contexts"),
    ]
    for x, name, color, note in cards:
        ax.add_patch(FancyBboxPatch((x, 1.1), 2.8, 2.6, boxstyle="round,pad=0.08,rounding_size=0.25",
                                    facecolor="white", edgecolor=color, lw=2))
        ax.text(x + 1.4, 3.1, name, ha="center", fontsize=14, fontweight="bold", color=color)
        ax.text(x + 1.4, 1.9, note, ha="center", va="center", fontsize=8.5, color="#475569")
    ax.text(6, 0.35, "Teleosts may carry a fourth paralog from the fish genome duplication.",
            ha="center", fontsize=8.5, color="#64748b")
    _save(fig, path)


def plot_homology_map(path: Path) -> None:
    """Sequence homology vs shared gap-junction job — keep the two ideas separate."""
    fig, ax = plt.subplots(figsize=(8.8, 4.4))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 7.4)
    ax.axis("off")
    ax.text(5, 7.0, "Sequence relationship", ha="center", fontsize=13, fontweight="bold", color="#1e3a4c")

    boxes = [
        (0.4, 4.3, "Innexins", "invertebrate\ngap junctions", "#0f766e"),
        (3.75, 4.3, "Pannexins", "chordate\nchannels (NGS)", "#c2410c"),
        (7.1, 4.3, "Connexins", "vertebrate\ngap junctions", "#0369a1"),
    ]
    for x, y, title, sub, color in boxes:
        ax.add_patch(
            FancyBboxPatch(
                (x, y),
                2.5,
                2.0,
                boxstyle="round,pad=0.06,rounding_size=0.2",
                facecolor="white",
                edgecolor=color,
                lw=2,
            )
        )
        ax.text(x + 1.25, y + 1.35, title, ha="center", fontsize=12, fontweight="bold", color=color)
        ax.text(x + 1.25, y + 0.55, sub, ha="center", fontsize=8.5, color="#475569")

    # Link innexin ↔ pannexin; label sits clearly below the boxes
    ax.annotate(
        "",
        xy=(3.75, 5.3),
        xytext=(2.9, 5.3),
        arrowprops=dict(arrowstyle="<->", color="#c2410c", lw=2.2),
    )
    ax.text(
        3.3,
        3.85,
        "homologous (same lineage)",
        ha="center",
        va="center",
        fontsize=9,
        color="#c2410c",
        fontweight="bold",
        bbox=dict(boxstyle="round,pad=0.25", facecolor="#fff7ed", edgecolor="#fdba74", linewidth=1),
    )

    ax.text(
        8.35,
        3.85,
        "no sequence homology\nto innexins / pannexins",
        ha="center",
        va="center",
        fontsize=8,
        color="#0369a1",
    )

    ax.add_patch(
        FancyBboxPatch(
            (0.5, 1.45),
            9.0,
            1.55,
            boxstyle="round,pad=0.04,rounding_size=0.15",
            facecolor="#f8fafc",
            edgecolor="#cbd5e1",
            lw=1,
        )
    )
    ax.text(
        5,
        2.55,
        "Same job, different genes",
        ha="center",
        fontsize=10,
        fontweight="bold",
        color="#334155",
    )
    ax.text(
        5,
        1.9,
        "Innexins and connexins both form gap junctions — that is analogy / niche, not homology.\n"
        "Pannexins are the retained innexin branch; they stay single-membrane channels (NGS).",
        ha="center",
        fontsize=8.2,
        color="#475569",
    )
    ax.text(
        5,
        0.5,
        "Innexin ↔ pannexin = homologous (lit.) · connexins = separate detectable sequence family in MMseqs panel",
        ha="center",
        fontsize=8.5,
        color="#1e3a4c",
    )
    _save(fig, path)


def write_html() -> None:
    html = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>A reading path for pannexins</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,560;9..144,700&family=Sora:wght@400;500;600&display=swap" rel="stylesheet">
<style>
:root {
  --ink:#1e2936; --muted:#5b6b78; --deep:#1e3a4c; --sea:#0e7490; --ember:#c2410c;
  --foam:#f3f6f8; --card:#fff; --line:rgba(30,41,54,.1);
}
* { box-sizing:border-box }
html { scroll-behavior:smooth }
body {
  margin:0; color:var(--ink); font-family:Sora,sans-serif; line-height:1.65;
  background:
    radial-gradient(900px 420px at 10% -10%, rgba(14,116,144,.14), transparent 55%),
    radial-gradient(800px 380px at 90% 0%, rgba(194,65,12,.10), transparent 50%),
    linear-gradient(180deg,#e7eef2 0%, var(--foam) 42%, #eef2f4 100%);
}
a { color:var(--sea) }
.hero {
  min-height:78vh; display:grid; align-items:end; padding:clamp(1.4rem,4vw,3rem); color:#f4fbff;
  background:linear-gradient(118deg, rgba(30,58,76,.95) 0%, rgba(14,116,144,.72) 52%, rgba(194,65,12,.55) 100%);
}
.brand { font-family:Fraunces,serif; font-size:clamp(2.2rem,6vw,4.2rem); font-weight:700; letter-spacing:-.03em; line-height:.95; margin:0 0 1rem }
.hero h1 { margin:0 0 .7rem; font-size:clamp(1.05rem,2vw,1.25rem); font-weight:500; max-width:40ch }
.hero .lede { margin:0 0 1.4rem; max-width:52ch; opacity:.93 }
.cta a {
  display:inline-block; text-decoration:none; font-weight:600; padding:.75rem 1.1rem; border-radius:999px;
  margin:.2rem .35rem .2rem 0; background:#f4fbff; color:var(--deep);
}
.cta a.ghost { background:transparent; color:#f4fbff; border:1px solid rgba(244,251,255,.45) }
.path {
  position:sticky; top:0; z-index:8; background:rgba(243,246,248,.92); backdrop-filter:blur(10px);
  border-bottom:1px solid var(--line); padding:.65rem 1rem; overflow:auto;
}
.path ol { display:flex; gap:.55rem; list-style:none; margin:0; padding:0; width:min(980px,100%); margin-inline:auto }
.path a { text-decoration:none; color:var(--muted); font-size:.78rem; white-space:nowrap; font-weight:600 }
.path a:hover { color:var(--sea) }
.wrap { width:min(980px, calc(100% - 2rem)); margin:0 auto; padding:2rem 0 4rem }
.chapter { padding:2.2rem 0; border-bottom:1px solid var(--line) }
.eyebrow { display:inline-block; font-size:.72rem; font-weight:700; letter-spacing:.08em; text-transform:uppercase; color:var(--ember); margin-bottom:.35rem }
.chapter h2 { font-family:Fraunces,serif; font-size:clamp(1.45rem,3vw,2rem); margin:0 0 .55rem; letter-spacing:-.02em }
.chapter .lead { color:var(--muted); max-width:62ch; margin:0 0 1.1rem }
.facts {
  display:grid; grid-template-columns:repeat(auto-fit,minmax(140px,1fr)); gap:.7rem; margin:1rem 0 1.3rem;
}
.facts div {
  background:var(--card); border:1px solid var(--line); border-radius:12px; padding:.85rem .9rem;
}
.facts b { display:block; font-family:Fraunces,serif; font-size:1.35rem; color:var(--deep) }
.facts span { color:var(--muted); font-size:.78rem }
.fig {
  margin:1rem 0; background:var(--card); border:1px solid var(--line); border-radius:16px; overflow:hidden;
}
.fig img { width:100%; display:block; cursor:zoom-in }
.fig figcaption { padding:.7rem .95rem; color:var(--muted); font-size:.82rem }
.card {
  display:block; text-decoration:none; color:inherit; margin-top:1rem;
  background:var(--card); border:1px solid var(--line); border-radius:14px; padding:1rem 1.1rem;
  box-shadow:0 6px 18px rgba(30,41,54,.05);
}
.card .kicker { font-size:.72rem; font-weight:700; letter-spacing:.06em; text-transform:uppercase; color:var(--ember); margin-bottom:.25rem }
.card h3 { font-family:Fraunces,serif; font-size:1.15rem; margin:0 0 .35rem; color:var(--deep) }
.card p { margin:0; color:var(--muted); font-size:.9rem; max-width:58ch }
.card .go { margin-top:.55rem; font-weight:600; color:var(--sea); font-size:.88rem }
.footer { margin-top:2rem; color:var(--muted); font-size:.85rem }
.footer a { margin-right:1rem }
.modal {
  position:fixed; inset:0; background:rgba(15,23,32,.82); display:none; align-items:center; justify-content:center; z-index:20; padding:1rem;
}
.modal.open { display:flex }
.modal img { max-width:min(1100px,96vw); max-height:90vh; border-radius:8px }
</style>
</head>
<body>
<header class="hero">
  <div>
    <p class="brand">Pannexins</p>
    <h1>The chordate side of the innexin story</h1>
    <p class="lede">
      The third family beside innexins and connexins: vertebrate-retained innexin homologs that work as
      single-membrane channels — docking blocked by constitutive N-linked glycosylation on the extracellular loops.
    </p>
    <div class="cta">
      <a href="#ch1">Start reading</a>
      <a class="ghost" href="tools.html">Tools &amp; plan</a>
    </div>
  </div>
</header>

<nav class="path" aria-label="Chapter path">
  <ol>
    <li><a href="#ch1">1 · What they are</a></li>
    <li><a href="#ch2">2 · NGS blockade</a></li>
    <li><a href="#ch3">3 · Chordate bottleneck</a></li>
    <li><a href="#ch4">4 · PANX1–3</a></li>
    <li><a href="#ch5">5 · Across clades</a></li>
    <li><a href="#ch6">6 · Homology</a></li>
    <li><a href="#ch7">7 · Results</a></li>
  </ol>
</nav>

<main class="wrap">
  <section class="chapter" id="ch1">
    <span class="eyebrow">Chapter 1</span>
    <h2>What pannexins are</h2>
    <p class="lead">
      The gap-junction path is mostly a binary invertebrate–vertebrate split (innexin × connexin).
      Pannexins are the third family: vertebrates kept ancestral innexin homologs, still with
      Cys motif 2 + 2 = 4 per monomer (not the connexin 3 + 3 = 6 pattern). Literature calls them
      heptameric single-membrane channels. The useful sequence marker here is a conserved
      extracellular N-X-S/T (Asn-X-Ser/Thr) site linked to non-docking behaviour — they do not
      form intercellular gap junctions.
    </p>
    <figure class="fig">
      <img src="figures/01_channel_vs_junction.png" alt="Pannexin channel versus classical gap junction" onclick="openModal(this.src)">
      <figcaption>Pannexin as a plasma-membrane channel vs a docking gap junction.</figcaption>
    </figure>
  </section>

  <section class="chapter" id="ch2">
    <span class="eyebrow">Chapter 2</span>
    <h2>The N-X-S/T sequence marker</h2>
    <p class="lead">
      Pannexins carry a conserved Asn-X-Ser/Thr (N-X-S/T) motif on their extracellular loops —
      the classical N-linked glycosylation (NGS) site pattern. That motif is a simple sequence
      filter separating non-docking pannexons from dockable innexin/connexin gap junctions.
      We use it as an annotation character, not as a sugar-chemistry argument.
    </p>
    <div class="facts">
      <div><b>N-X-S/T</b><span>EL motif</span></div>
      <div><b>Cys 2+2</b><span>innexin-like</span></div>
      <div><b>1 membrane</b><span>channel, not GJ</span></div>
      <div><b>≠ connexin</b><span>homolog of innexin</span></div>
    </div>
    <figure class="fig">
      <img src="figures/02_glycosylation_blockade.png" alt="N-X-S/T motif marks non-docking pannexins" onclick="openModal(this.src)">
      <figcaption>N-X-S/T on extracellular loops; innexin/connexin hemichannels remain dockable in the classical story.</figcaption>
    </figure>
  </section>

  <section class="chapter" id="ch3">
    <span class="eyebrow">Chapter 3</span>
    <h2>The early chordate bottleneck</h2>
    <p class="lead">
      Why keep a glycosylated innexin at all? Welzel &amp; Schuster 2022 (<em>eLife</em>) put the answer
      at the chordate origin: innexin diversity collapsed. Lancelets (<em>Branchiostoma</em>) have
      no connexins and only one conserved innexin with an extracellular NGS site. That remnant
      may not form functional gap junctions. The usual reading is that this bottleneck favoured
      later connexin diversification, while the glycosylated innexin lineage continued as today’s
      pannexins (PANX1–3 after vertebrate genome duplications).
    </p>
    <div class="facts">
      <div><b>1</b><span>innexin in lancelets (lit.)</span></div>
      <div><b>0</b><span>connexins in lancelets (lit.)</span></div>
      <div><b>NGS</b><span>on the remaining EL</span></div>
      <div><b>→ Cx</b><span>later GJ diversity</span></div>
    </div>
    <figure class="fig">
      <img src="figures/05_bottleneck_copies.png" alt="Early chordate vs mammal pannexin copy numbers" onclick="openModal(this.src)">
      <figcaption>
        Panel copy numbers for the literature bottleneck species versus sample mammals.
        Lamprey keeps one named protein here; <em>Branchiostoma</em> two (weak LOC labels);
        <em>Ciona</em> contributes no pannexin to the UniProt panel. The dated bottleneck timeline
        itself remains Welzel &amp; Schuster framing.
      </figcaption>
    </figure>
  </section>

  <section class="chapter" id="ch4">
    <span class="eyebrow">Chapter 4</span>
    <h2>Three named paralogs</h2>
    <p class="lead">
      Jawed vertebrates typically keep <strong>PANX1</strong>, <strong>PANX2</strong>, and <strong>PANX3</strong>
      (from early vertebrate genome duplications of that single glycosylated chordate innexin).
      PANX1 and PANX3 are closer; PANX2 is the outlier. Teleosts can carry an extra copy.
    </p>
    <div class="facts">
      <div><b>3</b><span>core paralogs</span></div>
      <div><b>PANX1</b><span>broad / ATP signalling</span></div>
      <div><b>PANX2</b><span>divergent / neural</span></div>
      <div><b>PANX3</b><span>closer to PANX1</span></div>
    </div>
    <figure class="fig">
      <img src="figures/03_three_paralogs.png" alt="PANX1 PANX2 PANX3 cards" onclick="openModal(this.src)">
      <figcaption>Working labels for the clade comparison that follows.</figcaption>
    </figure>
  </section>

  <section class="chapter" id="ch5">
    <span class="eyebrow">Chapter 5</span>
    <h2>Across clades: who keeps which paralogs</h2>
    <p class="lead">
      The curated UniProt panel is vertebrate-heavy on purpose — mammals, birds, reptiles, amphibians,
      ray- and lobe-finned fish, a cartilaginous fish — with early-chordate context from lamprey,
      lancelet, and tunicate. Across that panel, named PANX1 / PANX2 / PANX3 labels fill most jawed
      vertebrate rows as a triplet; early chordates stay sparse. PANX2 proteins run longer
      (median ~674 aa) than PANX1 / PANX3 (~426 / ~408 aa).
    </p>
    <div class="facts">
      <div><b>49</b><span>panel proteins</span></div>
      <div><b>17</b><span>species</span></div>
      <div><b>11</b><span>full PANX1–3 triplets</span></div>
      <div><b>9</b><span>clade bins</span></div>
    </div>
    <figure class="fig">
      <img src="figures/06_clade_type_heatmap.png" alt="Clade × PANX1/2/3 heatmap" onclick="openModal(this.src)">
      <figcaption>
        Clade × type counts from the curated panel. Mammals show balanced 1/2/3 columns;
        teleosts and cartilaginous fish are patchier; jawless / lancelet rows are thin.
      </figcaption>
    </figure>
    <figure class="fig">
      <img src="figures/07_type_mix.png" alt="PANX type mix stacked by clade" onclick="openModal(this.src)">
      <figcaption>Same data as stacked bars — which paralogs dominate each clade bin.</figcaption>
    </figure>
    <a class="card" href="../pannexin_clade_comparison/index.html">
      <div class="kicker">Full tables</div>
      <h3>Clade × type comparison</h3>
      <p>Length distributions, per-species copy numbers, and the annotation-status join — same shape as the connexin clade page.</p>
      <div class="go">Open clade comparison →</div>
    </a>
  </section>

  <section class="chapter" id="ch6">
    <span class="eyebrow">Chapter 6</span>
    <h2>Homology is the point</h2>
    <p class="lead">
      Unlike connexins, pannexins are homologous to innexins — the retained chordate branch after
      the diversity collapse. Connexins fill the vertebrate gap-junction niche as a separate sequence
      family. Full picture: innexin (invertebrate GJ), connexin (vertebrate GJ), pannexin
      (vertebrate channel, NGS-marked remnant).
    </p>
    <figure class="fig">
      <img src="figures/04_homology_map.png" alt="Innexin pannexin connexin relationship map" onclick="openModal(this.src)">
      <figcaption>Innexin ↔ pannexin homologous; connexins are a separate sequence family.</figcaption>
    </figure>
    <a class="card" href="../panx_vs_inx/index.html">
      <div class="kicker">Sequence check</div>
      <h3>Pannexin × innexin</h3>
      <p>Same MMseqs thresholds as innexin×connexin — here most pannexins hit at least one innexin.</p>
      <div class="go">Open homology contrast →</div>
    </a>
  </section>

  <section class="chapter" id="ch7">
    <span class="eyebrow">Chapter 7</span>
    <h2>What the panel shows</h2>
    <p class="lead">
      Taken together: a curated vertebrate/chordate UniProt set, clade × PANX1–3 occupancy,
      within-family similarity, a pannexin×innexin homology contrast, and a simple PANX tree
      with Ciona innexin outgroups. Weak LOC labels (lancelet, lamprey) still need the phylogeny
      page as a second look.
    </p>
    <div class="facts">
      <div><b>~90%</b><span>panx with ≥1 inx hit</span></div>
      <div><b>0</b><span>inx×cnx hits (parent)</span></div>
      <div><b>PANX1–3</b><span>jawed vertebrate core</span></div>
      <div><b>low</b><span>early-chordate copies</span></div>
    </div>
    <a class="card" href="../pannexin_insights/index.html">
      <div class="kicker">Gallery</div>
      <h3>Pannexin insights</h3>
      <p>Curated figures: homology, panx×inx, clades, rendered tree, exon notes.</p>
      <div class="go">Open insights →</div>
    </a>
    <a class="card" href="../pannexin_similarity/index.html">
      <div class="kicker">Within family</div>
      <h3>Pannexin similarity</h3>
      <p>MMseqs and composition space inside the curated panel.</p>
      <div class="go">Open similarity →</div>
    </a>
    <a class="card" href="../pannexin_phylogeny/index.html">
      <div class="kicker">Tree</div>
      <h3>PANX1 / PANX2 / PANX3 phylogeny</h3>
      <p>MAFFT + IQ-TREE on the panel; Ciona innexins as outgroups — no connexins.</p>
      <div class="go">Open phylogeny →</div>
    </a>
    <a class="card" href="../pannexin_annotation_status/index.html">
      <div class="kicker">Annotation</div>
      <h3>UniProt / NCBI status by clade</h3>
      <p>Same layered classifier used for innexins and connexins, for this species list.</p>
      <div class="go">Open annotation status →</div>
    </a>
  </section>

  <footer class="footer">
    <p>Pannexin reading path — chordate side of the innexin story.</p>
    <p>
      <a href="tools.html">Tools</a>
      <a href="../pannexin_clade_comparison/index.html">Clade × type</a>
      <a href="../panx_vs_inx/index.html">× innexin</a>
      <a href="../pannexin_phylogeny/index.html">Phylogeny</a>
      <a href="../../../../project/results/phylogenetic_story/index.html">Main phylogeny story</a>
      <a href="../../../../project/results/gap_junction_path/index.html">Gap-junction path</a>
      <a href="../../../../project/results/sources/index.html">Sources</a>
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
    html = r'''<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Tools &amp; plan — pannexins</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root { --ink:#1e2936; --muted:#5b6b78; --deep:#1e3a4c; --sea:#0e7490; --foam:#f3f6f8; --card:#fff; --line:rgba(30,41,54,.1); }
body { margin:0; font-family:Sora,sans-serif; color:var(--ink); background:linear-gradient(180deg,#e7eef2,var(--foam)); line-height:1.6 }
.hero { padding:1.6rem 1.4rem; color:#f4fbff; background:linear-gradient(118deg,#1e3a4c,#0e7490 60%,#c2410c) }
.hero p { font-family:Fraunces,serif; font-size:2rem; margin:0 0 .35rem }
.wrap { width:min(720px,calc(100% - 2rem)); margin:0 auto; padding:1.4rem 0 3rem }
h2 { font-family:Fraunces,serif; font-size:1.15rem }
li { background:var(--card); border:1px solid var(--line); border-radius:12px; padding:.7rem .9rem; margin:0 0 .45rem; list-style:none }
ul { padding:0 } li span { display:block; color:var(--muted); font-size:.82rem; margin-top:.15rem }
a { color:var(--sea) } .back { color:#f4fbff; font-weight:600; text-decoration:none }
</style>
</head>
<body>
<header class="hero">
  <p>Tools &amp; plan</p>
  <span>Starter stack for the pannexin project.</span>
  <div><a class="back" href="index.html">← path</a></div>
</header>
<main class="wrap">
  <h2>Intended toolkit (same family as the gap-junction thesis)</h2>
  <ul>
    <li><b>UniProt download</b><span>Reference PANX FASTAs per species.</span></li>
    <li><b>MMseqs2</b><span>Within-pannexin and pannexin→innexin (1 thread by default).</span></li>
    <li><b>MAFFT / IQ-TREE</b><span>Tree already on disk; light rebuilds reuse it.</span></li>
    <li><b>Python · Matplotlib · NumPy</b><span>Figures and composition embeddings.</span></li>
    <li><b>Neighbourhood searches</b><span>Optional later work — not required for the UniProt panel story.</span></li>
  </ul>
  <h2>Fixed points from the literature</h2>
  <ul>
    <li><b>Third family</b><span>Beside invertebrate innexins and vertebrate connexins: pannexins are retained innexin homologs in chordates/vertebrates.</span></li>
    <li><b>Heptamer</b><span>Panx1 cryo-EM: homo-heptamer (7) — not the innexon octamer (8) or connexon hexamer (6).</span></li>
    <li><b>Job</b><span>Single-membrane channels; marked by extracellular N-X-S/T; they do not form intercellular gap junctions.</span></li>
    <li><b>N-X-S/T</b><span>Conserved Asn-X-Ser/Thr motif on extracellular loops — sequence marker for the non-docking lineage.</span></li>
    <li><b>Early chordate bottleneck</b><span>Welzel &amp; Schuster 2022 (<em>eLife</em>): lancelets keep one N-X-S/T-marked innexin and no connexins.</span></li>
    <li><b>Cys motif</b><span>Same as innexins: 2 Cys per EL → 4 per monomer (not connexin 3+3=6).</span></li>
  </ul>
  <h2>Parameters</h2>
  <ul>
    <li><b>Not locked yet</b><span>When discovery and similarity builders land, we will pin the same style of explicit settings used on the gap-junction tools page (MMseqs <code>-s</code>/<code>-e</code>, cluster thresholds, IQ-TREE <code>-m MFP -bb 1000</code>, etc.).</span></li>
  </ul>
  <p><a href="index.html">← reading path</a></p>
</main>
</body>
</html>
'''
    (OUT / "tools.html").write_text(html, encoding="utf-8")


def _copy_clade_figures() -> None:
    """Embed clade-comparison plots into the path figures/ folder."""
    src_dir = RESULTS_DIR / "pannexin_clade_comparison" / "figures"
    mapping = {
        "07_bottleneck_focus.png": "05_bottleneck_copies.png",
        "02_clade_type_heatmap.png": "06_clade_type_heatmap.png",
        "03_type_mix_stacked.png": "07_type_mix.png",
    }
    for src_name, dest_name in mapping.items():
        src = src_dir / src_name
        if not src.exists():
            print(f"WARN missing {src} — run pannexin_clade_comparison_builder first")
            continue
        shutil.copy2(src, FIGS / dest_name)


def main() -> None:
    FIGS.mkdir(parents=True, exist_ok=True)
    plot_channel_vs_junction(FIGS / "01_channel_vs_junction.png")
    plot_glycosylation_blockade(FIGS / "02_glycosylation_blockade.png")
    plot_three_paralogs(FIGS / "03_three_paralogs.png")
    plot_homology_map(FIGS / "04_homology_map.png")
    _copy_clade_figures()
    write_html()
    write_tools_html()
    print(f"Done → {OUT / 'index.html'}")
    print(f"Tools → {OUT / 'tools.html'}")


if __name__ == "__main__":
    main()
