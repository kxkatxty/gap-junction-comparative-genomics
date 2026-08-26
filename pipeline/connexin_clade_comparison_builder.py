#!/usr/bin/env python3
"""Clade × connexin-type comparison (reference DB + discovery)."""

from __future__ import annotations

import html
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from pipeline.common import RESULTS_DIR, write_csv
from pipeline.connexin_common import (
    CLADE_COLORS,
    SUBFAMILY_COLORS,
    SUBFAMILY_FROM_TYPE,
    TYPE_COLORS,
    apply_best_hits,
    load_candidate_table_rows,
    load_discovery_fastas,
    load_exon_rows,
    load_reference_rows,
    summarize,
    write_fasta,
)

OUT = RESULTS_DIR / "connexin_clade_comparison"
FIGS = OUT / "figures"
WORK = OUT / "work"
MMSEQS = shutil.which("mmseqs") or str(
    Path.home() / "miniconda3/envs/synvoy_env/bin/mmseqs"
)

plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "DejaVu Sans",
    }
)


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def assign_discovery_types(refs: list[dict], disc: list[dict]) -> None:
    if not refs or not disc:
        return
    WORK.mkdir(parents=True, exist_ok=True)
    q = WORK / "discovery_query.fasta"
    t = WORK / "reference_targets.fasta"
    write_fasta(q, disc)
    write_fasta(t, refs)
    m8 = OUT / "discovery_vs_reference.m8"
    tmp = WORK / "search_tmp"
    tmp.mkdir(exist_ok=True)
    subprocess.run(
        [
            MMSEQS,
            "easy-search",
            str(q),
            str(t),
            str(m8),
            str(tmp),
            "--max-seqs",
            "5",
            "-e",
            "1e-3",
            "-v",
            "1",
            "--format-output",
            "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
        ],
        check=True,
    )
    apply_best_hits(disc, m8)


def plot_loci_by_clade(rows: list[dict], out: Path) -> None:
    sources = ["reference_db", "discovery"]
    labels = ["Reference DB", "Discovery"]
    colors = ["#9A3412", "#F59E0B"]
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    fig, ax = plt.subplots(figsize=(10, 5.2))
    x = np.arange(len(clades))
    bottom = np.zeros(len(clades))
    for src, lab, col in zip(sources, labels, colors):
        vals = np.array([sum(1 for r in rows if r["clade"] == c and r["source"] == src) for c in clades], dtype=float)
        ax.bar(x, vals, bottom=bottom, label=lab, color=col, width=0.72)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=25, ha="right")
    ax.set_ylabel("Connexin loci")
    ax.set_title("Connexin loci by vertebrate clade")
    ax.legend(frameon=False)
    _save(fig, out)


def plot_heatmap_clade_type(rows: list[dict], out: Path) -> None:
    types = [t for t, _ in Counter(r["reference_type"] for r in rows).most_common()]
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    mat = np.zeros((len(clades), len(types)))
    idx_c = {c: i for i, c in enumerate(clades)}
    idx_t = {t: i for i, t in enumerate(types)}
    for r in rows:
        mat[idx_c[r["clade"]], idx_t[r["reference_type"]]] += 1
    fig, ax = plt.subplots(figsize=(max(10, 0.65 * len(types) + 4), max(4.8, 0.5 * len(clades) + 2)))
    im = ax.imshow(mat, aspect="auto", cmap="YlOrBr")
    ax.set_xticks(range(len(types)))
    ax.set_xticklabels(types, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(clades)))
    ax.set_yticklabels(clades, fontsize=9)
    vmax = mat.max() or 1
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = int(mat[i, j])
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=7, color="#1c1917" if v < vmax * 0.65 else "white")
    ax.set_title("Clade × connexin type (gene name or best hit)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Locus count")
    _save(fig, out)


def plot_type_mix(rows: list[dict], out: Path) -> None:
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    types = [t for t, _ in Counter(r["reference_type"] for r in rows).most_common()]
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(clades))
    bottom = np.zeros(len(clades))
    for t in types:
        vals = np.array([sum(1 for r in rows if r["clade"] == c and r["reference_type"] == t) for c in clades], dtype=float)
        if vals.sum() == 0:
            continue
        ax.bar(x, vals, bottom=bottom, label=t, color=TYPE_COLORS.get(t, "#94A3B8"), width=0.75)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=25, ha="right")
    ax.set_ylabel("Loci")
    ax.set_title("Connexin-type mix within each vertebrate clade")
    ax.legend(ncol=3, fontsize=7, frameon=False, loc="upper right")
    _save(fig, out)


def plot_subfamily_by_clade(rows: list[dict], out: Path) -> None:
    sfs = ["alpha_GJA", "beta_GJB", "gamma_GJC", "delta_GJD", "epsilon_GJE", "unassigned"]
    labels = ["α GJA", "β GJB", "γ GJC", "δ GJD", "ε GJE", "unassigned"]
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    fig, ax = plt.subplots(figsize=(10, 5.2))
    x = np.arange(len(clades))
    bottom = np.zeros(len(clades))
    for sf, lab in zip(sfs, labels):
        vals = np.array([sum(1 for r in rows if r["clade"] == c and r["subfamily"] == sf) for c in clades], dtype=float)
        ax.bar(x, vals, bottom=bottom, label=lab, color=SUBFAMILY_COLORS.get(sf, "#94A3B8"), width=0.72)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=25, ha="right")
    ax.set_ylabel("Loci")
    ax.set_title("α / β / γ / δ / ε connexin subfamilies across clades")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, out)


def plot_classical_across_clades(rows: list[dict], out: Path) -> None:
    focus = ["GJA1", "GJA4", "GJA5", "GJB1", "GJB2", "GJB3", "GJC", "GJD"]
    clades = ["Mammal", "Bird", "Reptile", "Amphibian", "Ray-finned fish", "Cartilaginous fish"]
    clades = [c for c in clades if any(r["clade"] == c for r in rows)]
    fig, ax = plt.subplots(figsize=(11.5, 5.4))
    x = np.arange(len(focus))
    width = 0.13
    for i, clade in enumerate(clades):
        vals = [sum(1 for r in rows if r["reference_type"] == t and r["clade"] == clade) for t in focus]
        ax.bar(
            x + (i - len(clades) / 2) * width,
            vals,
            width=width,
            label=clade,
            color=CLADE_COLORS.get(clade, "#64748B"),
            alpha=0.92,
        )
    ax.set_xticks(x)
    ax.set_xticklabels(focus, rotation=20, ha="right")
    ax.set_ylabel("Locus count")
    ax.set_title("Classical connexin labels across vertebrate clades")
    ax.legend(fontsize=7, ncol=2, frameon=False)
    _save(fig, out)


def _boxplot(rows: list[dict], key: str, value: str, out: Path, title: str, ylabel: str, min_n: int = 2) -> None:
    groups: dict[str, list[float]] = defaultdict(list)
    for r in rows:
        v = r.get(value)
        if value == "reference_identity":
            n = float(v or 0)
            if n <= 0:
                continue
            groups[str(r[key])].append(n * (100 if n <= 1.5 else 1))
        else:
            n = float(v or 0)
            if n > 0:
                groups[str(r[key])].append(n)
    labels = [k for k, vals in sorted(groups.items(), key=lambda kv: -len(kv[1])) if len(vals) >= min_n]
    if not labels:
        return
    fig, ax = plt.subplots(figsize=(max(8, 0.55 * len(labels) + 3), 5))
    bp = ax.boxplot([groups[k] for k in labels], tick_labels=labels, patch_artist=True, showfliers=False)
    for patch, lab in zip(bp["boxes"], labels):
        patch.set_facecolor(CLADE_COLORS.get(lab) or TYPE_COLORS.get(lab) or SUBFAMILY_COLORS.get(lab) or "#64748B")
        patch.set_alpha(0.8)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
    _save(fig, out)


def plot_species_richness(rows: list[dict], out: Path) -> None:
    counts = Counter((r["clade"], r["organism"]) for r in rows if r["source"] == "discovery")
    items = sorted(counts.items(), key=lambda kv: -kv[1])[:25]
    if not items:
        return
    fig, ax = plt.subplots(figsize=(10, 6.8))
    y = np.arange(len(items))
    ax.barh(y, [c for _, c in items], color=[CLADE_COLORS.get(cl, "#64748B") for (cl, _), _ in items], alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([org for (_, org), _ in items], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("High-confidence discovery loci")
    ax.set_title("Top discovery species by connexin locus count")
    handles = [Patch(color=CLADE_COLORS[c], label=c) for c in sorted({cl for (cl, _), _ in items}) if c in CLADE_COLORS]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="lower right")
    _save(fig, out)


def plot_identity_vs_length(rows: list[dict], out: Path) -> None:
    disc = [r for r in rows if r["source"] == "discovery" and float(r.get("reference_identity") or 0) > 0 and r["protein_length"] > 0]
    if len(disc) < 5:
        return
    fig, ax = plt.subplots(figsize=(8, 5.4))
    for clade in sorted({r["clade"] for r in disc}):
        pts = [r for r in disc if r["clade"] == clade]
        ident = [float(r["reference_identity"]) * (100 if float(r["reference_identity"]) <= 1.5 else 1) for r in pts]
        ax.scatter(
            [r["protein_length"] for r in pts],
            ident,
            s=36,
            alpha=0.75,
            label=clade,
            c=CLADE_COLORS.get(clade, "#64748B"),
            edgecolors="white",
            linewidths=0.4,
        )
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("Best-hit identity (%)")
    ax.set_title("Discovery connexins: length vs identity by clade")
    ax.legend(fontsize=7, frameon=False)
    _save(fig, out)


def plot_exon_by_clade(exon_rows: list[dict], out: Path) -> None:
    groups: dict[str, list[float]] = defaultdict(list)
    for r in exon_rows:
        if r["mean_exon_count"] > 0:
            groups[r["clade"]].append(r["mean_exon_count"])
    labels = [k for k, v in sorted(groups.items(), key=lambda kv: -len(kv[1])) if len(v) >= 2]
    if not labels:
        return
    fig, ax = plt.subplots(figsize=(8.5, 5))
    bp = ax.boxplot([groups[k] for k in labels], tick_labels=labels, patch_artist=True, showfliers=False)
    for patch, lab in zip(bp["boxes"], labels):
        patch.set_facecolor(CLADE_COLORS.get(lab, "#64748B"))
        patch.set_alpha(0.8)
    ax.set_ylabel("Mean exons per gene")
    ax.set_title("Reference connexin exon architecture by clade")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=25, ha="right")
    _save(fig, out)


def plot_exon_by_type(exon_rows: list[dict], out: Path) -> None:
    groups: dict[str, list[float]] = defaultdict(list)
    for r in exon_rows:
        if r["mean_exon_count"] > 0:
            groups[r["subfamily"]].append(r["mean_exon_count"])
    labels = [k for k, v in sorted(groups.items(), key=lambda kv: -len(kv[1])) if v]
    if not labels:
        return
    fig, ax = plt.subplots(figsize=(8, 5))
    bp = ax.boxplot([groups[k] for k in labels], tick_labels=labels, patch_artist=True, showfliers=False)
    for patch, lab in zip(bp["boxes"], labels):
        patch.set_facecolor(SUBFAMILY_COLORS.get(lab, "#64748B"))
        patch.set_alpha(0.8)
    ax.set_ylabel("Mean exons")
    ax.set_title("Exon count by connexin subfamily (α/β/γ/δ)\n(δ GJD spikes often = GFF fusion with flanks, not true 22-exon genes)")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=20, ha="right")
    ax.text(
        0.02,
        0.98,
        "Connexin ORFs are typically intron-free (~2 exons with 5′ UTR).\n"
        "Extreme exon counts often fuse compact connexins with intron-rich flanks\n"
        "(e.g. DLGAP3 / SMIM12 near GJA4).",
        transform=ax.transAxes,
        va="top",
        fontsize=7.5,
        color="#64748B",
    )
    _save(fig, out)


def write_tables(rows: list[dict]) -> None:
    slim = [{k: v for k, v in r.items() if k != "seq"} for r in rows]
    write_csv(OUT / "all_connexin_loci_annotated.csv", slim)
    clade_rows = []
    for clade in sorted({r["clade"] for r in rows}):
        sub = [r for r in rows if r["clade"] == clade]
        ids = [
            float(r["reference_identity"]) * (100 if float(r["reference_identity"]) <= 1.5 else 1)
            for r in sub
            if float(r.get("reference_identity") or 0) > 0 and r["source"] != "reference_db"
        ]
        lens = [r["protein_length"] for r in sub if r["protein_length"] > 0]
        clade_rows.append(
            {
                "clade": clade,
                "n_loci": len(sub),
                "n_species": len({r["species_dir"] for r in sub}),
                "n_reference": sum(1 for r in sub if r["source"] == "reference_db"),
                "n_discovery": sum(1 for r in sub if r["source"] == "discovery"),
                "median_length_aa": f"{np.median(lens):.0f}" if lens else "",
                "median_identity_pct": f"{np.median(ids):.1f}" if ids else "",
                "top_types": "; ".join(f"{k}:{v}" for k, v in Counter(r["reference_type"] for r in sub).most_common(5)),
            }
        )
    write_csv(OUT / "clade_summary.csv", clade_rows)
    type_rows = []
    for t in sorted({r["reference_type"] for r in rows}):
        sub = [r for r in rows if r["reference_type"] == t]
        type_rows.append(
            {
                "reference_type": t,
                "subfamily": SUBFAMILY_FROM_TYPE.get(t, "unassigned"),
                "n_loci": len(sub),
                "n_clades": len({r["clade"] for r in sub}),
                "clades": "; ".join(sorted({r["clade"] for r in sub})),
                "n_discovery": sum(1 for r in sub if r["source"] == "discovery"),
                "n_reference": sum(1 for r in sub if r["source"] == "reference_db"),
            }
        )
    write_csv(OUT / "reference_type_summary.csv", type_rows)


def key_findings(rows: list[dict]) -> list[str]:
    findings = []
    refs = [r for r in rows if r["source"] == "reference_db"]
    disc = [r for r in rows if r["source"] == "discovery"]
    findings.append(
        f"Dataset: {len(refs)} reference connexins + {len(disc)} high-confidence discovery loci "
        f"across {len({r['species_dir'] for r in rows})} species."
    )
    if disc:
        top_t, n_t = Counter(r["reference_type"] for r in disc).most_common(1)[0]
        findings.append(
            f"Discovery recoveries most often match {top_t} ({n_t}/{len(disc)}). "
            "Labels are best-hit similarity, not proven 1:1 orthology."
        )
    # alpha vs beta presence across clades
    clades = sorted({r["clade"] for r in rows})
    both = [
        c
        for c in clades
        if any(r["subfamily"] == "alpha_GJA" and r["clade"] == c for r in rows)
        and any(r["subfamily"] == "beta_GJB" and r["clade"] == c for r in rows)
    ]
    if both:
        findings.append(
            f"α (GJA) and β (GJB) subfamilies co-occur in: {', '.join(both)} — "
            "the dual-locus architecture (GJA1 vs β-cluster) is a vertebrate-wide pattern, not mammal-only."
        )
    if disc:
        top_c, n_c = Counter(r["clade"] for r in disc).most_common(1)[0]
        findings.append(f"Among discovery loci, {top_c} contributes the most ({n_c}).")

    def med_len(clade: str) -> float | None:
        vals = [r["protein_length"] for r in rows if r["clade"] == clade and r["protein_length"] > 80]
        return float(np.median(vals)) if vals else None

    for a, b in (("Mammal", "Amphibian"), ("Bird", "Ray-finned fish"), ("Amphibian", "Cartilaginous fish")):
        ma, mb = med_len(a), med_len(b)
        if ma and mb:
            findings.append(f"Median protein length: {a} {ma:.0f} aa vs {b} {mb:.0f} aa.")
    findings.append(
        "Reference set covers α/β/γ/δ/ε classes (GJA, GJB, GJC, GJD, GJE). "
        "GJA4 sits in the conserved β-neighborhood; GJA1 is a separate locus."
    )
    findings.append(
        "Note: δ GJD exon-count outliers (up to ~22) are often automated GFF fusions "
        "with intron-rich flanks (e.g. DLGAP3/SMIM12 near GJA4), not real 22-exon connexins."
    )
    return findings


def build_html(rows: list[dict], findings: list[str], fig_plan: list[tuple[str, str]]) -> None:
    stats = summarize(rows)
    finding_html = "".join(f"<li>{html.escape(f)}</li>" for f in findings)
    cards = []
    for fname, caption in fig_plan:
        if not (FIGS / fname).exists():
            continue
        cards.append(
            f"""<article class="card">
  <img src="figures/{html.escape(fname)}" alt="{html.escape(caption)}" loading="lazy" onclick="openModal(this.src)">
  <h3>{html.escape(caption)}</h3>
</article>"""
        )
    clade_table = ["<tr><th>Clade</th><th>Loci</th><th>Species</th><th>Ref</th><th>Discovery</th></tr>"]
    for clade, n in stats["by_clade"].most_common():
        sub = [r for r in rows if r["clade"] == clade]
        clade_table.append(
            "<tr>"
            f"<td>{html.escape(clade)}</td><td>{n}</td>"
            f"<td>{len({r['species_dir'] for r in sub})}</td>"
            f"<td>{sum(1 for r in sub if r['source']=='reference_db')}</td>"
            f"<td>{sum(1 for r in sub if r['source']=='discovery')}</td>"
            "</tr>"
        )
    type_table = ["<tr><th>Type</th><th>Subfamily</th><th>Loci</th><th>Clades</th></tr>"]
    for t, n in stats["by_type"].most_common():
        type_table.append(
            f"<tr><td>{html.escape(t)}</td><td>{html.escape(SUBFAMILY_FROM_TYPE.get(t,'unassigned'))}</td>"
            f"<td>{n}</td><td>{len({r['clade'] for r in rows if r['reference_type']==t})}</td></tr>"
        )
    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Connexin Clade × Type Comparison</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root {{ --ink:#2a1408; --muted:#7a5340; --deep:#7c2d12; --amber:#c2410c; --bg:#fbf4ee; }}
*{{box-sizing:border-box}}
body{{margin:0;font-family:Sora,sans-serif;color:var(--ink);background:
 radial-gradient(900px 420px at 0% 0%, rgba(234,88,12,.14), transparent 55%),
 linear-gradient(180deg,#f8ebe3,var(--bg));line-height:1.55}}
.hero{{min-height:78vh;display:grid;align-items:end;padding:clamp(1.4rem,4vw,3rem);color:#fffaf6;
 background:linear-gradient(120deg,rgba(124,45,18,.94),rgba(194,65,12,.78) 55%,rgba(217,119,6,.55))}}
.brand{{font-family:Fraunces,serif;font-size:clamp(2.4rem,6.5vw,4.6rem);margin:0 0 .7rem;line-height:.95}}
.hero h1{{font-size:clamp(1.05rem,2vw,1.3rem);font-weight:500;max-width:42ch;margin:0 0 .7rem}}
.hero p{{max-width:52ch;opacity:.92}}
.cta a{{display:inline-block;text-decoration:none;background:#fffaf6;color:var(--deep);padding:.8rem 1.15rem;border-radius:999px;font-weight:700;margin:.4rem .4rem 0 0}}
.cta a.ghost{{background:transparent;color:#fffaf6;border:1px solid rgba(255,250,246,.45)}}
.wrap{{width:min(1160px,calc(100% - 2rem));margin:0 auto;padding:2rem 0 4rem}}
.stats{{display:grid;grid-template-columns:repeat(4,1fr);gap:.8rem;margin-top:-2rem}}
.stat{{background:#fff;border-radius:16px;padding:1rem;box-shadow:0 10px 28px rgba(42,20,8,.06)}}
.stat b{{display:block;font-family:Fraunces,serif;font-size:1.8rem;color:var(--deep)}}
.stat span{{color:var(--muted);font-size:.84rem}}
h2{{font-family:Fraunces,serif}}
.findings{{background:#fff;border-radius:16px;padding:1.1rem 1.3rem;border-left:4px solid #d97706}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:1rem}}
.card{{background:#fff;border-radius:16px;padding:.85rem;box-shadow:0 8px 22px rgba(42,20,8,.05)}}
.card img{{width:100%;border-radius:10px;cursor:zoom-in}}
.card h3{{font-size:.95rem;margin:.7rem 0 0}}
table{{width:100%;border-collapse:collapse;background:#fff;border-radius:12px;overflow:hidden;font-size:.88rem}}
th,td{{padding:.55rem .7rem;border-bottom:1px solid #f3e0d4;text-align:left}}
th{{background:#7c2d12;color:#fffaf6}}
.modal{{display:none;position:fixed;inset:0;background:rgba(0,0,0,.88);z-index:20;align-items:center;justify-content:center;padding:1.5rem;cursor:zoom-out}}
.modal.open{{display:flex}} .modal img{{max-width:95vw;max-height:90vh;border-radius:8px}}
@media(max-width:800px){{.stats{{grid-template-columns:1fr 1fr}}}}
</style></head>
<body>
<header class="hero"><div>
<p class="brand">Clade × Connexin</p>
<h1>Compare vertebrate clades and classical connexin labels in one place</h1>
<p>Reference proteins plus high-confidence discovery loci — sliced by clade, gene type (GJA1, GJB2, GJC…) and α/β/γ/δ/ε subfamily.</p>
<div class="cta">
<a href="#figures">See diagrams</a>
<a class="ghost" href="../innexin_clade_comparison/index.html">Innexin counterpart</a>
</div>
</div></header>
<main class="wrap">
<div class="stats">
<div class="stat"><b>{stats['n_loci']}</b><span>annotated connexin loci</span></div>
<div class="stat"><b>{stats['n_species']}</b><span>species represented</span></div>
<div class="stat"><b>{len(stats['by_clade'])}</b><span>clade bins</span></div>
<div class="stat"><b>{len(stats['by_type'])}</b><span>type labels</span></div>
</div>
<section><h2>Key results</h2>
<div class="findings"><ol>{finding_html}</ol></div>
<p style="color:#7a5340;max-width:70ch">A discovery locus labelled GJA1-like means its <em>best reference hit</em> was GJA1. Orthology still needs synteny (GJA1 vs the β-cluster around GJA4).</p>
</section>
<section><h2>Clade overview</h2><table>{''.join(clade_table)}</table></section>
<section><h2>Connexin labels</h2><table>{''.join(type_table)}</table></section>
<section id="figures"><h2>Diagrams</h2>
<div class="grid">{''.join(cards)}</div>
</section>
<section>
<p>
<a href="all_connexin_loci_annotated.csv">all loci</a> ·
<a href="clade_summary.csv">clade summary</a> ·
<a href="reference_type_summary.csv">type summary</a> ·
<a href="../connexin_similarity/index.html">MMseqs + embeddings</a> ·
<a href="../connexin_insights/index.html">Curated insights</a> ·
<a href="../family_comparison/index.html">vs innexin</a>
</p>
</section>
</main>
<div id="modal" class="modal" onclick="this.classList.remove('open')"><img id="modal-img" alt=""></div>
<script>
function openModal(src){{document.getElementById('modal-img').src=src;document.getElementById('modal').classList.add('open')}}
document.addEventListener('keydown',e=>{{if(e.key==='Escape')document.getElementById('modal').classList.remove('open')}});
</script>
</body></html>"""
    (OUT / "index.html").write_text(page, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    refs = load_reference_rows()
    table = load_candidate_table_rows()
    disc_seq = load_discovery_fastas(high_conf_only=True)
    print(f"Refs={len(refs)} table_kept={len(table)} discovery_fasta_hc={len(disc_seq)}")
    assign_discovery_types(refs, disc_seq)

    # Merge: references + discovery sequences (typed) + table rows without seq (for counts if missing fasta)
    have = {(r["species_dir"], r["candidate_id"]) for r in disc_seq}
    extra = []
    for r in table:
        if (r["species_dir"], r["candidate_id"]) in have:
            continue
        if "high_confidence" not in (r.get("rank_category") or "").lower() and "possible_connexin" not in (r.get("rank_category") or "").lower():
            continue
        extra.append(r)
    rows = refs + disc_seq + extra
    write_tables(rows)
    exon_rows = load_exon_rows()

    plot_loci_by_clade(rows, FIGS / "01_loci_by_clade_source.png")
    plot_heatmap_clade_type(rows, FIGS / "02_heatmap_clade_x_type.png")
    plot_type_mix(rows, FIGS / "03_type_mix_per_clade.png")
    plot_classical_across_clades(rows, FIGS / "04_classical_across_clades.png")
    plot_subfamily_by_clade(rows, FIGS / "05_subfamily_by_clade.png")
    _boxplot(rows, "clade", "protein_length", FIGS / "06_length_by_clade.png", "Protein length by clade", "Length (aa)")
    _boxplot(rows, "reference_type", "protein_length", FIGS / "07_length_by_type.png", "Protein length by connexin type", "Length (aa)")
    _boxplot(
        [r for r in rows if r["source"] == "discovery"],
        "clade",
        "reference_identity",
        FIGS / "08_identity_by_clade.png",
        "Best-hit identity by clade (discovery)",
        "Identity (%)",
    )
    _boxplot(
        [r for r in rows if r["source"] == "discovery"],
        "reference_type",
        "reference_identity",
        FIGS / "09_identity_by_type.png",
        "Best-hit identity by connexin type",
        "Identity (%)",
    )
    plot_species_richness(rows, FIGS / "10_species_richness.png")
    plot_identity_vs_length(rows, FIGS / "11_identity_vs_length.png")
    plot_exon_by_clade(exon_rows, FIGS / "13_exon_by_clade.png")
    plot_exon_by_type(exon_rows, FIGS / "14_exon_by_subfamily.png")

    fig_plan = [
        ("01_loci_by_clade_source.png", "Loci by clade, stacked by source"),
        ("02_heatmap_clade_x_type.png", "Heatmap: clade × connexin type"),
        ("03_type_mix_per_clade.png", "Type mix within clades"),
        ("04_classical_across_clades.png", "GJA1/GJA4/GJB… across clades"),
        ("05_subfamily_by_clade.png", "α/β/γ/δ/ε subfamilies by clade"),
        ("06_length_by_clade.png", "Protein length by clade"),
        ("07_length_by_type.png", "Protein length by connexin type"),
        ("08_identity_by_clade.png", "Best-hit identity by clade"),
        ("09_identity_by_type.png", "Best-hit identity by type"),
        ("10_species_richness.png", "Top discovery species"),
        ("11_identity_vs_length.png", "Length vs identity scatter"),
        ("13_exon_by_clade.png", "Exon architecture by clade"),
        ("14_exon_by_subfamily.png", "Exons by α/β/γ/δ (δ spikes = annotation fusion)"),
    ]
    findings = key_findings(rows)
    (OUT / "key_findings.md").write_text("# Key findings\n\n" + "\n".join(f"- {f}" for f in findings) + "\n", encoding="utf-8")
    build_html(rows, findings, fig_plan)
    print(f"Wrote {OUT}/index.html  loci={len(rows)}")


if __name__ == "__main__":
    main()
