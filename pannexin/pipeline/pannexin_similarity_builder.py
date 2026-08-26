#!/usr/bin/env python3
"""Within-pannexin MMseqs similarity + type/clade overview.

Outputs → project/results/pannexin_similarity/
"""

from __future__ import annotations

import html
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pipeline.common import RESULTS_DIR, TYPE_COLORS, read_csv_rows, write_csv
from pipeline.limits import n_threads
from pipeline.pannexin_panel_builder import FASTA as PANEL_FASTA
from pipeline.pannexin_panel_builder import META as PANEL_META
from pipeline.pannexin_panel_builder import build as build_panel

OUT = RESULTS_DIR / "pannexin_similarity"
FIGS = OUT / "figures"
WORK = OUT / "work"
FASTA = OUT / "all_pannexins.fasta"
META = OUT / "sequence_metadata.csv"
M8 = OUT / "mmseqs_allvsall.m8"

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
    }
)


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _ensure_panel() -> list[dict]:
    build_panel()
    shutil.copy2(PANEL_FASTA, FASTA)
    shutil.copy2(PANEL_META, META)
    return read_csv_rows(META)


def _run_mmseqs(fasta: Path, m8: Path) -> None:
    if not Path(MMSEQS).exists():
        raise FileNotFoundError(f"mmseqs not found: {MMSEQS}")
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    db = WORK / "db"
    res = WORK / "res"
    tmp = WORK / "tmp"
    tmp.mkdir()
    subprocess.run([MMSEQS, "createdb", str(fasta), str(db)], check=True, capture_output=True)
    subprocess.run(
        [
            MMSEQS,
            "search",
            str(db),
            str(db),
            str(res),
            str(tmp),
            "-e",
            "1e-3",
            "--threads",
            str(n_threads(1)),
            "-s",
            "5.0",
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            MMSEQS,
            "convertalis",
            str(db),
            str(db),
            str(res),
            str(m8),
            "--format-output",
            "query,target,pident,alnlen,evalue,bits",
        ],
        check=True,
        capture_output=True,
    )


def _parse_m8(path: Path) -> list[dict]:
    rows = []
    if not path.exists():
        return rows
    for line in path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        q, t, pident, alnlen, evalue, bits = parts[:6]
        if q == t:
            continue
        rows.append(
            {
                "query": q,
                "target": t,
                "pident": float(pident),
                "alnlen": int(float(alnlen)),
                "evalue": float(evalue),
                "bits": float(bits),
            }
        )
    return rows


def _type_of(seq_id: str, meta: dict[str, dict]) -> str:
    return meta.get(seq_id, {}).get("panx_type", "other/unknown")


def plot_length_by_type(rows: list[dict], path: Path) -> None:
    groups: dict[str, list[int]] = defaultdict(list)
    for r in rows:
        groups[r["panx_type"]].append(int(r["length"]))
    labels = [t for t in ("PANX1", "PANX2", "PANX3", "other/unknown") if t in groups]
    data = [groups[t] for t in labels]
    fig, ax = plt.subplots(figsize=(7.2, 4.0))
    bp = ax.boxplot(data, tick_labels=labels, patch_artist=True)
    for patch, lab in zip(bp["boxes"], labels):
        patch.set_facecolor(TYPE_COLORS.get(lab, "#94A3B8"))
        patch.set_alpha(0.75)
    ax.set_ylabel("protein length (aa)")
    ax.set_title("Pannexin length by paralog")
    _save(fig, path)


def plot_identity_heatmap(hits: list[dict], meta: dict[str, dict], path: Path) -> None:
    # Median best reciprocal identity within/between types
    types = ["PANX1", "PANX2", "PANX3", "other/unknown"]
    present = [t for t in types if any(_type_of(h["query"], meta) == t for h in hits) or any(
        r.get("panx_type") == t for r in meta.values()
    )]
    best: dict[tuple[str, str], list[float]] = defaultdict(list)
    for h in hits:
        if h["alnlen"] < 80:
            continue
        a = _type_of(h["query"], meta)
        b = _type_of(h["target"], meta)
        best[(a, b)].append(h["pident"])
    mat = np.full((len(present), len(present)), np.nan)
    for i, a in enumerate(present):
        for j, b in enumerate(present):
            vals = best.get((a, b), [])
            if vals:
                mat[i, j] = float(np.median(vals))
    fig, ax = plt.subplots(figsize=(5.6, 4.8))
    im = ax.imshow(mat, cmap="YlOrBr", vmin=20, vmax=90)
    ax.set_xticks(range(len(present)))
    ax.set_yticks(range(len(present)))
    ax.set_xticklabels(present, rotation=30, ha="right")
    ax.set_yticklabels(present)
    for i in range(len(present)):
        for j in range(len(present)):
            v = mat[i, j]
            if not np.isnan(v):
                ax.text(j, i, f"{v:.0f}", ha="center", va="center", fontsize=9)
    ax.set_title("Median MMseqs %id between pannexin types")
    fig.colorbar(im, ax=ax, fraction=0.046, label="% identity")
    _save(fig, path)


def plot_clade_counts(rows: list[dict], path: Path) -> None:
    counts: dict[str, int] = defaultdict(int)
    for r in rows:
        counts[r["clade"]] += 1
    labels = sorted(counts, key=lambda k: -counts[k])
    fig, ax = plt.subplots(figsize=(7.5, 3.8))
    ax.bar(labels, [counts[k] for k in labels], color="#C2410C", alpha=0.85)
    ax.set_ylabel("proteins")
    ax.set_title("Pannexin panel by clade")
    ax.tick_params(axis="x", rotation=25)
    _save(fig, path)


def write_html(rows: list[dict], n_hits: int) -> None:
    n = len(rows)
    n_sp = len({r["organism"] for r in rows})
    by_type = defaultdict(int)
    for r in rows:
        by_type[r["panx_type"]] += 1
    type_bits = "".join(
        f"<li><strong>{html.escape(t)}</strong>: {c}</li>" for t, c in sorted(by_type.items())
    )
    page = f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pannexin similarity</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root {{ --ink:#12202a; --muted:#4d6570; --ember:#c2410c; --line:rgba(18,32,42,.12); }}
body {{ margin:0; font-family:Sora,sans-serif; color:var(--ink); background:#f7f3ef;
  background-image:radial-gradient(900px 420px at 90% 0%, rgba(194,65,12,.12), transparent 55%); }}
.wrap {{ width:min(860px, calc(100% - 2rem)); margin:0 auto; padding:2rem 0 4rem }}
h1 {{ font-family:Fraunces,serif; font-size:clamp(1.6rem,4vw,2.2rem); margin:0 0 .4rem }}
.lede {{ color:var(--muted); max-width:58ch }}
.fig {{ background:#fff; border:1px solid var(--line); border-radius:14px; padding:.8rem; margin:1rem 0 }}
.fig img {{ width:100%; display:block; border-radius:8px }}
.cap {{ font-size:.88rem; color:var(--muted); margin:.5rem 0 0 }}
.footer a {{ margin-right:.85rem; color:var(--ember) }}
</style></head><body>
<main class="wrap">
  <h1>Pannexin similarity</h1>
  <p class="lede">Within-family MMseqs2 for the curated chordate/vertebrate panel
  (<strong>{n}</strong> proteins, <strong>{n_sp}</strong> species; {n_hits} off-diagonal hits at e≤1e-3).</p>
  <ul>{type_bits}</ul>
  <div class="fig"><img src="figures/01_length_by_type.png" alt="length by type">
    <p class="cap">PANX2 is typically longer (extended C-terminus); PANX1/PANX3 cluster near ~400 aa.</p></div>
  <div class="fig"><img src="figures/02_type_identity.png" alt="type identity">
    <p class="cap">Median pairwise % identity within and between paralog classes (aln ≥80 aa).</p></div>
  <div class="fig"><img src="figures/03_clade_counts.png" alt="clade counts">
    <p class="cap">Panel composition by clade (UniProt download; false hits filtered).</p></div>
  <p class="footer">
    <a href="../pannexin_path/index.html">Pannexin path</a>
    <a href="../panx_vs_inx/index.html">Pannexin × innexin</a>
    <a href="../pannexin_panel/sequence_metadata.csv">Panel table</a>
  </p>
</main></body></html>
"""
    (OUT / "index.html").write_text(page, encoding="utf-8")


def build() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    rows = _ensure_panel()
    meta = {r["seq_id"]: r for r in rows}
    _run_mmseqs(FASTA, M8)
    hits = _parse_m8(M8)
    write_csv(OUT / "mmseqs_hits_summary.csv", [
        {
            "n_proteins": str(len(rows)),
            "n_offdiag_hits": str(len(hits)),
            "median_pident": f"{np.median([h['pident'] for h in hits]):.1f}" if hits else "",
        }
    ])
    plot_length_by_type(rows, FIGS / "01_length_by_type.png")
    plot_identity_heatmap(hits, meta, FIGS / "02_type_identity.png")
    plot_clade_counts(rows, FIGS / "03_clade_counts.png")
    write_html(rows, len(hits))
    print(f"Done → {OUT / 'index.html'} ({len(rows)} proteins, {len(hits)} hits)")


if __name__ == "__main__":
    build()
