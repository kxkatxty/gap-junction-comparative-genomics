#!/usr/bin/env python3
"""Clade × PANX1/2/3 comparison for the curated pannexin panel.

Mirrors innexin/connexin clade_comparison pages (tables + figures + HTML)
using UniProt panel proteins only — no SynVoy, no discovery genomes.

Outputs → project/results/pannexin_clade_comparison/
"""

from __future__ import annotations

import html
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from pipeline.common import (
    CLADE_COLORS,
    RESULTS_DIR,
    TYPE_COLORS,
    read_csv_rows,
    write_csv,
)

OUT = RESULTS_DIR / "pannexin_clade_comparison"
FIGS = OUT / "figures"
PANEL_META = RESULTS_DIR / "pannexin_panel" / "sequence_metadata.csv"
STATUS_CSV = RESULTS_DIR / "pannexin_annotation_status" / "species_status.csv"

TYPE_ORDER = ("PANX1", "PANX2", "PANX3", "other/unknown")
BOTTLENECK_CLADES = {"Jawless vertebrate", "Tunicate / lancelet"}

plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
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


def load_rows() -> list[dict]:
    rows = []
    for r in read_csv_rows(PANEL_META):
        try:
            length = int(float(r.get("length") or 0))
        except (TypeError, ValueError):
            length = 0
        try:
            cys = int(float(r.get("cys_count") or 0))
        except (TypeError, ValueError):
            cys = 0
        clade = (r.get("clade") or "Other").strip() or "Other"
        panx_type = (r.get("panx_type") or "other/unknown").strip() or "other/unknown"
        rows.append(
            {
                "seq_id": r.get("seq_id", ""),
                "organism": r.get("organism", ""),
                "organism_key": r.get("organism_key", ""),
                "gene": r.get("gene", ""),
                "accession": r.get("accession", ""),
                "panx_type": panx_type,
                "clade": clade,
                "length": length,
                "cys_count": cys,
                "source": "reference_panel",
                "bottleneck_clade": "yes" if clade in BOTTLENECK_CLADES else "no",
            }
        )
    return rows


def write_tables(rows: list[dict]) -> None:
    write_csv(
        OUT / "all_pannexin_loci_annotated.csv",
        [
            {
                **{k: str(v) for k, v in r.items()},
            }
            for r in rows
        ],
    )

    clade_rows = []
    for clade in sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c)):
        sub = [r for r in rows if r["clade"] == clade]
        lens = [r["length"] for r in sub if r["length"] > 0]
        type_counts = Counter(r["panx_type"] for r in sub)
        clade_rows.append(
            {
                "clade": clade,
                "n_loci": str(len(sub)),
                "n_species": str(len({r["organism"] for r in sub})),
                "median_length_aa": f"{float(np.median(lens)):.0f}" if lens else "",
                "n_panx1": str(type_counts.get("PANX1", 0)),
                "n_panx2": str(type_counts.get("PANX2", 0)),
                "n_panx3": str(type_counts.get("PANX3", 0)),
                "n_other": str(type_counts.get("other/unknown", 0)),
                "species": "; ".join(sorted({r["organism"] for r in sub})),
                "bottleneck": "yes" if clade in BOTTLENECK_CLADES else "no",
            }
        )
    write_csv(OUT / "clade_summary.csv", clade_rows)

    type_rows = []
    for t in TYPE_ORDER:
        sub = [r for r in rows if r["panx_type"] == t]
        if not sub:
            continue
        lens = [r["length"] for r in sub if r["length"] > 0]
        type_rows.append(
            {
                "panx_type": t,
                "n_loci": str(len(sub)),
                "n_clades": str(len({r["clade"] for r in sub})),
                "n_species": str(len({r["organism"] for r in sub})),
                "median_length_aa": f"{float(np.median(lens)):.0f}" if lens else "",
                "clades": "; ".join(sorted({r["clade"] for r in sub})),
            }
        )
    write_csv(OUT / "reference_type_summary.csv", type_rows)

    # species × type copy-number matrix (bottleneck automation)
    species_rows = []
    by_sp: dict[str, list[dict]] = defaultdict(list)
    for r in rows:
        by_sp[r["organism"]].append(r)
    for org, items in sorted(by_sp.items(), key=lambda kv: (kv[1][0]["clade"], kv[0])):
        tc = Counter(r["panx_type"] for r in items)
        species_rows.append(
            {
                "organism": org,
                "clade": items[0]["clade"],
                "n_loci": str(len(items)),
                "n_panx1": str(tc.get("PANX1", 0)),
                "n_panx2": str(tc.get("PANX2", 0)),
                "n_panx3": str(tc.get("PANX3", 0)),
                "n_other": str(tc.get("other/unknown", 0)),
                "has_full_triplet": "yes"
                if tc.get("PANX1", 0) and tc.get("PANX2", 0) and tc.get("PANX3", 0)
                else "no",
                "bottleneck_clade": items[0]["bottleneck_clade"],
            }
        )
    write_csv(OUT / "species_copy_number.csv", species_rows)

    # clade × type count matrix
    matrix_rows = []
    clades = [r["clade"] for r in clade_rows]
    for clade in clades:
        row = {"clade": clade}
        for t in TYPE_ORDER:
            row[t] = str(sum(1 for r in rows if r["clade"] == clade and r["panx_type"] == t))
        matrix_rows.append(row)
    write_csv(OUT / "clade_by_type_matrix.csv", matrix_rows)


def plot_loci_by_clade(rows: list[dict], out: Path) -> None:
    counts = Counter(r["clade"] for r in rows)
    clades = [c for c, _ in counts.most_common()]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    x = np.arange(len(clades))
    vals = [counts[c] for c in clades]
    colors = [CLADE_COLORS.get(c, "#94A3B8") for c in clades]
    ax.bar(x, vals, color=colors, width=0.72)
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=28, ha="right")
    ax.set_ylabel("Proteins in panel")
    ax.set_title("Pannexin panel: loci by clade")
    _save(fig, out)


def plot_clade_type_heatmap(rows: list[dict], out: Path) -> None:
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    types = [t for t in TYPE_ORDER if any(r["panx_type"] == t for r in rows)]
    mat = np.array(
        [[sum(1 for r in rows if r["clade"] == c and r["panx_type"] == t) for t in types] for c in clades],
        dtype=float,
    )
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    im = ax.imshow(mat, cmap="YlOrBr", aspect="auto")
    ax.set_xticks(np.arange(len(types)))
    ax.set_yticks(np.arange(len(clades)))
    ax.set_xticklabels(types)
    ax.set_yticklabels(clades)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            ax.text(j, i, int(mat[i, j]), ha="center", va="center", color="#1e293b", fontsize=9)
    ax.set_title("Clade × PANX type counts")
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04, label="proteins")
    _save(fig, out)


def plot_type_mix(rows: list[dict], out: Path) -> None:
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    fig, ax = plt.subplots(figsize=(10, 5.2))
    x = np.arange(len(clades))
    bottom = np.zeros(len(clades))
    for t in TYPE_ORDER:
        vals = np.array([sum(1 for r in rows if r["clade"] == c and r["panx_type"] == t) for c in clades], dtype=float)
        if vals.sum() == 0:
            continue
        ax.bar(x, vals, bottom=bottom, label=t, color=TYPE_COLORS.get(t, "#94A3B8"), width=0.72)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=28, ha="right")
    ax.set_ylabel("Proteins")
    ax.set_title("PANX1 / PANX2 / PANX3 mix within each clade")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, out)


def plot_length_by_type(rows: list[dict], out: Path) -> None:
    types = [t for t in TYPE_ORDER if any(r["panx_type"] == t and r["length"] > 0 for r in rows)]
    data = [[r["length"] for r in rows if r["panx_type"] == t and r["length"] > 0] for t in types]
    fig, ax = plt.subplots(figsize=(7.5, 4.8))
    bp = ax.boxplot(data, patch_artist=True, widths=0.55)
    ax.set_xticklabels(types)
    for patch, t in zip(bp["boxes"], types):
        patch.set_facecolor(TYPE_COLORS.get(t, "#94A3B8"))
        patch.set_alpha(0.75)
    ax.set_ylabel("Length (aa)")
    ax.set_title("Protein length by PANX type")
    _save(fig, out)


def plot_length_by_clade(rows: list[dict], out: Path) -> None:
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    data = [[r["length"] for r in rows if r["clade"] == c and r["length"] > 0] for c in clades]
    fig, ax = plt.subplots(figsize=(10, 5))
    bp = ax.boxplot(data, patch_artist=True, widths=0.55)
    for patch, c in zip(bp["boxes"], clades):
        patch.set_facecolor(CLADE_COLORS.get(c, "#94A3B8"))
        patch.set_alpha(0.75)
    ax.set_xticklabels(clades, rotation=28, ha="right")
    ax.set_ylabel("Length (aa)")
    ax.set_title("Protein length by clade")
    _save(fig, out)


def plot_species_richness(rows: list[dict], out: Path) -> None:
    counts = Counter(r["organism"] for r in rows)
    items = sorted(counts.items(), key=lambda kv: -kv[1])
    org_clade = {r["organism"]: r["clade"] for r in rows}
    fig, ax = plt.subplots(figsize=(9, 7))
    y = np.arange(len(items))
    vals = [c for _, c in items]
    colors = [CLADE_COLORS.get(org_clade.get(o, "Other"), "#94A3B8") for o, _ in items]
    ax.barh(y, vals, color=colors, alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels([o for o, _ in items], fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Panel proteins")
    ax.set_title("Copy number per species (UniProt panel)")
    handles = [
        Patch(color=CLADE_COLORS[c], label=c)
        for c in sorted({org_clade[o] for o, _ in items})
        if c in CLADE_COLORS
    ]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="lower right")
    _save(fig, out)


def plot_bottleneck_focus(rows: list[dict], out: Path) -> None:
    """Automated lancelet / lamprey / tunicate copy-number check."""
    focus = [r for r in rows if r["bottleneck_clade"] == "yes"]
    # also show a typical jawed vertebrate for scale
    mammal = [r for r in rows if r["clade"] == "Mammal"]
    orgs_focus = sorted({r["organism"] for r in focus})
    ref_orgs = sorted({r["organism"] for r in mammal})[:3]
    orgs = orgs_focus + [o for o in ref_orgs if o not in orgs_focus]
    if not orgs:
        return
    types = list(TYPE_ORDER)
    fig, ax = plt.subplots(figsize=(9, 4.8))
    x = np.arange(len(orgs))
    width = 0.18
    for i, t in enumerate(types):
        vals = [sum(1 for r in rows if r["organism"] == o and r["panx_type"] == t) for o in orgs]
        ax.bar(x + (i - 1.5) * width, vals, width=width, label=t, color=TYPE_COLORS.get(t, "#94A3B8"))
    ax.set_xticks(x)
    ax.set_xticklabels(orgs, rotation=20, ha="right", fontsize=8)
    ax.set_ylabel("Panel proteins")
    ax.set_title("Bottleneck check: early chordates vs sample mammals")
    ax.legend(frameon=False, fontsize=8, ncol=4)
    _save(fig, out)


def plot_cys_by_type(rows: list[dict], out: Path) -> None:
    types = [t for t in TYPE_ORDER if any(r["panx_type"] == t for r in rows)]
    data = [[r["cys_count"] for r in rows if r["panx_type"] == t] for t in types]
    fig, ax = plt.subplots(figsize=(7, 4.5))
    bp = ax.boxplot(data, patch_artist=True, widths=0.55)
    ax.set_xticklabels(types)
    for patch, t in zip(bp["boxes"], types):
        patch.set_facecolor(TYPE_COLORS.get(t, "#94A3B8"))
        patch.set_alpha(0.75)
    ax.set_ylabel("Cys count")
    ax.set_title("Cysteine count by PANX type (panel)")
    _save(fig, out)


def key_findings(rows: list[dict]) -> list[str]:
    findings = []
    n_sp = len({r["organism"] for r in rows})
    findings.append(
        f"Curated UniProt panel: {len(rows)} pannexin proteins in {n_sp} species across "
        f"{len({r['clade'] for r in rows})} clade bins."
    )
    tc = Counter(r["panx_type"] for r in rows)
    findings.append(
        "Type mix: "
        + ", ".join(f"{t}={tc.get(t, 0)}" for t in TYPE_ORDER if tc.get(t, 0))
        + "."
    )
    # triplet completeness
    by_sp: dict[str, Counter] = defaultdict(Counter)
    for r in rows:
        by_sp[r["organism"]][r["panx_type"]] += 1
    full = sum(1 for c in by_sp.values() if c["PANX1"] and c["PANX2"] and c["PANX3"])
    findings.append(
        f"{full}/{len(by_sp)} panel species carry at least one named PANX1, PANX2 and PANX3."
    )
    # bottleneck
    for org in sorted({r["organism"] for r in rows if r["bottleneck_clade"] == "yes"}):
        sub = [r for r in rows if r["organism"] == org]
        findings.append(
            f"Bottleneck species {org}: {len(sub)} panel protein(s) "
            f"({', '.join(f'{k}:{v}' for k, v in Counter(r['panx_type'] for r in sub).most_common())})."
        )
    # length
    for t in ("PANX1", "PANX2", "PANX3"):
        lens = [r["length"] for r in rows if r["panx_type"] == t and r["length"] > 0]
        if lens:
            findings.append(f"Median length {t}: {float(np.median(lens)):.0f} aa (n={len(lens)}).")
    # Ciona absence
    status = read_csv_rows(STATUS_CSV)
    for r in status:
        if r.get("organism") == "Ciona intestinalis":
            findings.append(
                f"Ciona intestinalis: UniProt/NCBI status={r.get('category_short', '?')}, "
                f"panel proteins={r.get('n_panel_proteins', '0')} "
                "(expect innexin, not pannexin)."
            )
    findings.append(
        "Orthology caveat: LOC / weak UniProt labels (lancelet, lamprey) are best-effort types — "
        "confirm with the pannexin phylogeny page."
    )
    return findings


def build_html(rows: list[dict], findings: list[str], fig_names: list[tuple[str, str]]) -> None:
    finding_html = "".join(f"<li>{html.escape(f)}</li>" for f in findings)
    cards = []
    for fname, caption in fig_names:
        if not (FIGS / fname).exists():
            continue
        cards.append(
            f"""
<article class="card">
  <img src="figures/{html.escape(fname)}" alt="{html.escape(caption)}" loading="lazy" onclick="openModal(this.src)">
  <h3>{html.escape(caption)}</h3>
</article>"""
        )

    clade_table = ["<tr><th>Clade</th><th>Loci</th><th>Species</th><th>PANX1</th><th>PANX2</th><th>PANX3</th><th>other</th><th>Median aa</th></tr>"]
    for r in read_csv_rows(OUT / "clade_summary.csv"):
        clade_table.append(
            "<tr>"
            f"<td>{html.escape(r['clade'])}</td>"
            f"<td>{html.escape(r['n_loci'])}</td>"
            f"<td>{html.escape(r['n_species'])}</td>"
            f"<td>{html.escape(r['n_panx1'])}</td>"
            f"<td>{html.escape(r['n_panx2'])}</td>"
            f"<td>{html.escape(r['n_panx3'])}</td>"
            f"<td>{html.escape(r['n_other'])}</td>"
            f"<td>{html.escape(r['median_length_aa'])}</td>"
            "</tr>"
        )

    copy_table = [
        "<tr><th>Species</th><th>Clade</th><th>n</th><th>1/2/3</th><th>full triplet?</th><th>bottleneck?</th></tr>"
    ]
    for r in read_csv_rows(OUT / "species_copy_number.csv"):
        copy_table.append(
            "<tr>"
            f"<td>{html.escape(r['organism'])}</td>"
            f"<td>{html.escape(r['clade'])}</td>"
            f"<td>{html.escape(r['n_loci'])}</td>"
            f"<td>{html.escape(r['n_panx1'])}/{html.escape(r['n_panx2'])}/{html.escape(r['n_panx3'])}</td>"
            f"<td>{html.escape(r['has_full_triplet'])}</td>"
            f"<td>{html.escape(r['bottleneck_clade'])}</td>"
            "</tr>"
        )

    n_loci = len(rows)
    n_species = len({r["organism"] for r in rows})
    n_clades = len({r["clade"] for r in rows})
    n_types = len({r["panx_type"] for r in rows})

    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pannexin Clade × Type Comparison</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600;9..144,700&family=Sora:wght@400;550;700&display=swap" rel="stylesheet">
<style>
:root {{
  --ink:#1c1917; --muted:#57534e; --ember:#c2410c; --deep:#7c2d12; --sand:#d4a017; --bg:#f7f3ef;
}}
*{{box-sizing:border-box}}
body{{margin:0;font-family:Sora,sans-serif;color:var(--ink);background:
 radial-gradient(900px 420px at 0% 0%, rgba(194,65,12,.10), transparent 55%),
 radial-gradient(800px 380px at 100% 10%, rgba(15,118,110,.10), transparent 50%),
 linear-gradient(180deg,#faf6f1,var(--bg));line-height:1.55}}
.hero{{min-height:72vh;display:grid;align-items:end;padding:clamp(1.4rem,4vw,3rem);color:#fff7ed;
 background:linear-gradient(120deg,rgba(124,45,18,.94),rgba(194,65,12,.72) 55%,rgba(15,118,110,.45)),
 radial-gradient(circle at 75% 25%,#ea580c,#7c2d12 60%)}}
.brand{{font-family:Fraunces,serif;font-size:clamp(2.2rem,6vw,4.2rem);margin:0 0 .7rem;line-height:.95;letter-spacing:-.03em}}
.hero h1{{font-size:clamp(1.05rem,2vw,1.3rem);font-weight:500;max-width:42ch;margin:0 0 .7rem}}
.hero p{{max-width:54ch;opacity:.92;margin:0 0 1.4rem}}
.cta a{{display:inline-block;text-decoration:none;background:#fff7ed;color:var(--deep);padding:.8rem 1.15rem;border-radius:999px;font-weight:700;margin-right:.5rem}}
.cta a.ghost{{background:transparent;color:#fff7ed;border:1px solid rgba(255,247,237,.45)}}
.wrap{{width:min(1160px,calc(100% - 2rem));margin:0 auto;padding:2rem 0 4rem}}
.stats{{display:grid;grid-template-columns:repeat(4,1fr);gap:.8rem;margin-top:-2rem;position:relative;z-index:2}}
.stat{{background:rgba(255,255,255,.92);border:1px solid rgba(28,25,23,.08);border-radius:16px;padding:1rem;box-shadow:0 10px 28px rgba(28,25,23,.06)}}
.stat b{{display:block;font-family:Fraunces,serif;font-size:1.8rem;color:var(--deep)}}
.stat span{{color:var(--muted);font-size:.84rem}}
section{{margin:2.4rem 0}}
h2{{font-family:Fraunces,serif;font-size:1.7rem;margin:0 0 .6rem}}
.lead{{color:var(--muted);max-width:70ch}}
.findings{{background:#fff;border-radius:16px;padding:1.1rem 1.3rem;border-left:4px solid var(--sand);box-shadow:0 8px 24px rgba(28,25,23,.05)}}
.findings ol{{margin:.4rem 0 0;padding-left:1.2rem}}
.findings li{{margin:.35rem 0;color:var(--muted)}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:1rem}}
.card{{background:rgba(255,255,255,.92);border:1px solid rgba(28,25,23,.08);border-radius:16px;padding:.85rem;box-shadow:0 8px 22px rgba(28,25,23,.05)}}
.card img{{width:100%;border-radius:10px;cursor:zoom-in;background:#f8fafc}}
.card h3{{font-size:.95rem;margin:.7rem 0 0;font-weight:600}}
table{{width:100%;border-collapse:collapse;background:#fff;border-radius:12px;overflow:hidden;box-shadow:0 6px 18px rgba(28,25,23,.05);font-size:.88rem}}
th,td{{padding:.55rem .7rem;border-bottom:1px solid #e7e5e4;text-align:left}}
th{{background:#7c2d12;color:#fff7ed;font-weight:600}}
.links a{{margin-right:1rem;color:var(--ember)}}
.modal{{display:none;position:fixed;inset:0;background:rgba(0,0,0,.88);z-index:20;align-items:center;justify-content:center;padding:1.5rem;cursor:zoom-out}}
.modal.open{{display:flex}}
.modal img{{max-width:95vw;max-height:90vh;border-radius:8px}}
@media(max-width:800px){{.stats{{grid-template-columns:1fr 1fr}}.hero{{min-height:62vh}}}}
</style>
</head>
<body>
<header class="hero">
  <div>
    <p class="brand">Clade × Pannexin</p>
    <h1>PANX1 / PANX2 / PANX3 across vertebrates and early chordates</h1>
    <p>Same shape as the connexin clade page: clade × type counts, copy number, lengths —
    built from the curated UniProt panel.</p>
    <div class="cta">
      <a href="#figures">See diagrams</a>
      <a class="ghost" href="../pannexin_annotation_status/index.html">Annotation status</a>
    </div>
  </div>
</header>
<main class="wrap">
  <div class="stats">
    <div class="stat"><b>{n_loci}</b><span>panel proteins</span></div>
    <div class="stat"><b>{n_species}</b><span>species</span></div>
    <div class="stat"><b>{n_clades}</b><span>clade bins</span></div>
    <div class="stat"><b>{n_types}</b><span>type labels</span></div>
  </div>

  <section>
    <h2>Key results</h2>
    <div class="findings"><ol>{finding_html}</ol></div>
    <p class="lead" style="margin-top:1rem">
      Types come from UniProt gene / description labels via the same PANX1–3 classifier used in the panel builder.
      Weak LOC names (lancelet, lamprey) need the phylogeny page for a second look.
    </p>
  </section>

  <section id="figures">
    <h2>Diagrams</h2>
    <div class="grid">{"".join(cards)}</div>
  </section>

  <section>
    <h2>Clade summary</h2>
    <table>{"".join(clade_table)}</table>
  </section>

  <section>
    <h2>Copy number (bottleneck automation)</h2>
    <p class="lead">Per-species PANX1/2/3 counts. Early-chordate rows are flagged for the literature bottleneck check.</p>
    <table>{"".join(copy_table)}</table>
  </section>

  <section class="links">
    <p>
      <a href="../pannexin_annotation_status/index.html">Annotation status</a>
      <a href="../pannexin_similarity/index.html">Similarity</a>
      <a href="../pannexin_phylogeny/index.html">Phylogeny</a>
      <a href="../panx_vs_inx/index.html">× innexin</a>
      <a href="../pannexin_path/index.html">Path</a>
      <a href="../pannexin_curator/index.html">Gaps</a>
    </p>
  </section>
</main>
<div class="modal" id="modal" onclick="this.classList.remove('open')"><img id="modal-img" alt=""></div>
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

    findings_md = OUT / "key_findings.md"
    findings_md.write_text("# Pannexin clade × type — key findings\n\n" + "\n".join(f"- {f}" for f in findings) + "\n", encoding="utf-8")


def build() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rows = load_rows()
    if not rows:
        raise SystemExit(f"No panel metadata at {PANEL_META} — run pannexin_panel_builder first.")

    write_tables(rows)

    fig_plan = [
        ("01_loci_by_clade.png", "Loci by clade", plot_loci_by_clade),
        ("02_clade_type_heatmap.png", "Clade × type heatmap", plot_clade_type_heatmap),
        ("03_type_mix_stacked.png", "Type mix within clades", plot_type_mix),
        ("04_length_by_type.png", "Length by PANX type", plot_length_by_type),
        ("05_length_by_clade.png", "Length by clade", plot_length_by_clade),
        ("06_species_copy_number.png", "Copy number per species", plot_species_richness),
        ("07_bottleneck_focus.png", "Lancelet / lamprey bottleneck check", plot_bottleneck_focus),
        ("08_cys_by_type.png", "Cys count by type", plot_cys_by_type),
    ]
    fig_names: list[tuple[str, str]] = []
    for fname, caption, fn in fig_plan:
        fn(rows, FIGS / fname)
        fig_names.append((fname, caption))

    findings = key_findings(rows)
    build_html(rows, findings, fig_names)
    print(f"Clade comparison: {len(rows)} proteins → {OUT / 'index.html'}")


if __name__ == "__main__":
    build()
