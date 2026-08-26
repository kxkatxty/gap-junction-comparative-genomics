#!/usr/bin/env python3
"""Pannexin × innexin contrast (homology expected; no connexin comparison).

Uses the curated pannexin panel FASTA and the parent-project innexin panel
(`../project/results/innexin_similarity/all_innexins.fasta`).

Outputs → project/results/panx_vs_inx/
"""

from __future__ import annotations

import html
import math
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pipeline.common import PROJECT_ROOT, RESULTS_DIR, write_csv
from pipeline.common import read_csv_rows
from pipeline.limits import n_threads
from pipeline.pannexin_panel_builder import FASTA as PANX_PANEL_FASTA
from pipeline.pannexin_panel_builder import META as PANX_META
from pipeline.pannexin_panel_builder import build as build_panel

OUT = RESULTS_DIR / "panx_vs_inx"
FIGS = OUT / "figures"
WORK = OUT / "work"

# Parent gap-junction project (one level above pannexin/)
PARENT_RESULTS = PROJECT_ROOT.parent / "project" / "results"
INX_FASTA = PARENT_RESULTS / "innexin_similarity" / "all_innexins.fasta"

PANX = "#C2410C"
INX = "#0F766E"
MUTED = "#64748B"

MMSEQS = shutil.which("mmseqs") or str(
    Path.home() / "miniconda3/envs/synvoy_env/bin/mmseqs"
)
EVALUE = "1e-3"
MIN_ALN = 80

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
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def read_fasta(path: Path, family: str) -> list[dict]:
    rows: list[dict] = []
    seq_id, chunks = "", []
    if not path.exists():
        return rows
    with path.open(encoding="utf-8") as handle:
        for line in handle:
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if seq_id:
                    rows.append(_feat(seq_id, "".join(chunks), family))
                seq_id = line[1:].split()[0]
                chunks = []
            else:
                chunks.append(line.upper())
        if seq_id:
            rows.append(_feat(seq_id, "".join(chunks), family))
    return rows


def _feat(seq_id: str, seq: str, family: str) -> dict:
    seq = "".join(a for a in seq if a.isalpha())
    n = max(len(seq), 1)
    return {
        "seq_id": seq_id,
        "family": family,
        "seq": seq,
        "length": len(seq),
        "cys_per_100": 100.0 * seq.count("C") / n,
        "gravy": sum(KD.get(a, 0.0) for a in seq) / n,
        "charged_frac": sum(seq.count(a) for a in "DEKR") / n,
    }


def kmer_svd(seqs: list[str], k: int = 3, n_comp: int = 6) -> np.ndarray:
    vocab_count: Counter[str] = Counter()
    for seq in seqs:
        for i in range(max(0, len(seq) - k + 1)):
            km = seq[i : i + k]
            if "X" in km:
                continue
            vocab_count[km] += 1
    vocab = [km for km, _ in vocab_count.most_common(600)]
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
    X -= X.mean(axis=0, keepdims=True)
    if X.shape[1] < 2 or X.shape[0] < 2:
        return np.zeros((len(seqs), 2))
    _, _, Vt = np.linalg.svd(X, full_matrices=False)
    comps = min(n_comp, Vt.shape[0], X.shape[1])
    return X @ Vt[:comps].T


def run_cross_mmseqs(panx: Path, inx: Path, m8: Path) -> list[dict]:
    if WORK.exists():
        shutil.rmtree(WORK)
    WORK.mkdir(parents=True)
    qdb, tdb, res, tmp = WORK / "qdb", WORK / "tdb", WORK / "res", WORK / "tmp"
    tmp.mkdir()
    subprocess.run([MMSEQS, "createdb", str(panx), str(qdb)], check=True, capture_output=True)
    subprocess.run([MMSEQS, "createdb", str(inx), str(tdb)], check=True, capture_output=True)
    subprocess.run(
        [
            MMSEQS, "search", str(qdb), str(tdb), str(res), str(tmp),
            "-e", EVALUE, "--threads", str(n_threads(1)), "-s", "5.7",
        ],
        check=True,
        capture_output=True,
    )
    subprocess.run(
        [
            MMSEQS, "convertalis", str(qdb), str(tdb), str(res), str(m8),
            "--format-output", "query,target,pident,alnlen,evalue,bits",
        ],
        check=True,
        capture_output=True,
    )
    hits = []
    for line in m8.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 6:
            continue
        q, t, pident, alnlen, evalue, bits = parts[:6]
        aln = int(float(alnlen))
        if aln < MIN_ALN:
            continue
        hits.append(
            {
                "query": q,
                "target": t,
                "pident": float(pident),
                "alnlen": aln,
                "evalue": float(evalue),
                "bits": float(bits),
            }
        )
    return hits


def plot_length(panx: list[dict], inx: list[dict], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    bp = ax.boxplot(
        [[r["length"] for r in panx], [r["length"] for r in inx]],
        tick_labels=["pannexin", "innexin"],
        patch_artist=True,
    )
    for patch, col in zip(bp["boxes"], [PANX, INX]):
        patch.set_facecolor(col)
        patch.set_alpha(0.75)
    ax.set_ylabel("length (aa)")
    ax.set_title("Protein length: pannexin vs innexin panels")
    _save(fig, path)


def plot_embedding(panx: list[dict], inx: list[dict], path: Path) -> None:
    all_rows = panx + inx
    coords = kmer_svd([r["seq"] for r in all_rows], k=3, n_comp=4)
    fig, ax = plt.subplots(figsize=(7.2, 5.2))
    n_p = len(panx)
    ax.scatter(coords[n_p:, 0], coords[n_p:, 1], s=18, c=INX, alpha=0.55, label=f"innexin (n={len(inx)})")
    ax.scatter(coords[:n_p, 0], coords[:n_p, 1], s=42, c=PANX, alpha=0.9, label=f"pannexin (n={len(panx)})", edgecolors="#7c2d12", linewidths=0.4)
    ax.set_xlabel("3-mer SVD 1")
    ax.set_ylabel("3-mer SVD 2")
    ax.set_title("Sequence-composition space (k=3)")
    ax.legend(frameon=False)
    _save(fig, path)


def plot_hit_identity(hits: list[dict], path: Path) -> None:
    fig, ax = plt.subplots(figsize=(7.0, 4.0))
    if hits:
        vals = [h["pident"] for h in hits]
        ax.hist(vals, bins=20, color=PANX, alpha=0.8, edgecolor="white")
        ax.axvline(np.median(vals), color=INX, lw=1.6, label=f"median {np.median(vals):.1f}%")
        ax.legend(frameon=False)
    else:
        ax.text(0.5, 0.5, "no cross-family hits", ha="center", transform=ax.transAxes)
    ax.set_xlabel("MMseqs % identity (aln ≥80 aa, e≤1e-3)")
    ax.set_ylabel("hit count")
    ax.set_title("Cross-family hits: pannexin query → innexin target")
    _save(fig, path)


def plot_coverage(panx: list[dict], hits: list[dict], path: Path) -> None:
    """Fraction of panx proteins with ≥1 innexin hit."""
    with_hit = {h["query"] for h in hits}
    n = len(panx)
    n_hit = sum(1 for r in panx if r["seq_id"] in with_hit)
    fig, ax = plt.subplots(figsize=(5.5, 3.6))
    ax.bar(["with ≥1 innexin hit", "no hit"], [n_hit, n - n_hit], color=[INX, MUTED], alpha=0.85)
    ax.set_ylabel("pannexin proteins")
    ax.set_title(f"Cross-family recovery ({n_hit}/{n} = {100*n_hit/max(n,1):.0f}%)")
    _save(fig, path)


def write_html(stats: dict) -> None:
    page = f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pannexin × innexin</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600;700&family=Sora:wght@400;500;600&display=swap" rel="stylesheet">
<style>
:root {{ --ink:#12202a; --muted:#4d6570; --panx:#c2410c; --inx:#0f766e; --line:rgba(18,32,42,.12); --foam:#f6f3ef; }}
body {{ margin:0; font-family:Sora,sans-serif; color:var(--ink); line-height:1.65;
  background:radial-gradient(900px 420px at 8% -5%, rgba(15,118,110,.14), transparent 55%),
             radial-gradient(800px 380px at 95% 0%, rgba(194,65,12,.12), transparent 50%),
             linear-gradient(180deg,#efe8e1 0%, var(--foam) 40%, #f3f0eb 100%); }}
.wrap {{ width:min(860px, calc(100% - 2rem)); margin:0 auto; padding:2rem 0 4rem }}
.brand {{ font-family:Fraunces,serif; font-size:clamp(1.7rem,4vw,2.4rem); margin:0 0 .35rem; letter-spacing:-.02em }}
.lede {{ color:var(--muted); max-width:58ch; margin:0 0 1.2rem }}
.facts {{ display:grid; grid-template-columns:repeat(4,1fr); gap:.65rem; margin:0 0 1.4rem }}
.fact {{ background:#fff; border:1px solid var(--line); border-radius:12px; padding:.75rem }}
.fact b {{ display:block; font-family:Fraunces,serif; font-size:1.35rem }}
.fact span {{ font-size:.75rem; color:var(--muted) }}
.take {{ border-left:3px solid var(--inx); background:rgba(15,118,110,.07); padding:.75rem .9rem; border-radius:0 12px 12px 0; margin:0 0 1.2rem }}
.fig {{ background:#fff; border:1px solid var(--line); border-radius:14px; padding:.85rem; margin:0 0 1rem }}
.fig img {{ width:100%; display:block; border-radius:8px }}
.cap {{ font-size:.88rem; color:var(--muted); margin:.55rem 0 0 }}
.closing {{ margin-top:1.8rem; padding:1.2rem 1.25rem; border-radius:14px; background:linear-gradient(135deg,#1e3a4c,var(--panx)); color:#fff }}
.footer {{ margin-top:1.4rem; color:var(--muted); font-size:.88rem }}
.footer a {{ margin-right:.85rem; color:var(--panx) }}
@media (max-width:700px) {{ .facts {{ grid-template-columns:1fr 1fr }} }}
</style></head><body>
<main class="wrap">
  <p class="brand">Pannexin × innexin</p>
  <p class="lede">Same comparison style as the innexin × connexin page, but only against innexins —
  the family pannexins are expected to resemble. Connexins are intentionally not in this contrast.</p>

  <div class="facts">
    <div class="fact"><b>{stats['n_panx']}</b><span>pannexin proteins</span></div>
    <div class="fact"><b>{stats['n_inx']}</b><span>innexin proteins</span></div>
    <div class="fact"><b>{stats['n_hits']}</b><span>cross hits (e≤1e-3, ≥80 aa)</span></div>
    <div class="fact"><b>{stats['pct_with_hit']}</b><span>panx with ≥1 inx hit</span></div>
  </div>

  <p class="take"><strong>Finding.</strong> {html.escape(stats['verdict'])}</p>

  <div class="fig"><img src="figures/01_lengths.png" alt="lengths">
    <p class="cap">Length distributions for the two panels.</p></div>
  <div class="fig"><img src="figures/02_kmer_space.png" alt="kmer space">
    <p class="cap">3-mer composition SVD. Overlap is expected if pannexins are innexin-related sequences;
    separation would argue against detectable relatedness in this feature space.</p></div>
  <div class="fig"><img src="figures/03_cross_identity.png" alt="cross identity">
    <p class="cap">MMseqs % identity for pannexin→innexin hits (aln ≥80 aa, e≤1e-3).</p></div>
  <div class="fig"><img src="figures/04_hit_coverage.png" alt="coverage">
    <p class="cap">How many pannexin proteins recover at least one innexin hit under the same thresholds
    used for innexin×connexin (where cross hits were essentially zero).</p></div>

  <section class="closing">
    <p>{html.escape(stats['closing'])}</p>
  </section>
  <p class="footer">
    <a href="../pannexin_path/index.html">Pannexin path</a>
    <a href="../pannexin_similarity/index.html">Within-pannexin similarity</a>
    <a href="../../../../project/results/inx_vs_cnx/index.html">Innexin × connexin (parent)</a>
  </p>
</main></body></html>
"""
    (OUT / "index.html").write_text(page, encoding="utf-8")


def build() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    build_panel()
    if not INX_FASTA.exists():
        raise FileNotFoundError(
            f"Innexin FASTA missing: {INX_FASTA}\n"
            "Build the parent innexin similarity panel first."
        )

    panx = read_fasta(PANX_PANEL_FASTA, "pannexin")
    inx = read_fasta(INX_FASTA, "innexin")
    m8 = OUT / "panx_vs_inx.m8"
    hits = run_cross_mmseqs(PANX_PANEL_FASTA, INX_FASTA, m8)
    write_csv(OUT / "cross_hits.csv", [
        {k: str(v) for k, v in h.items()} for h in hits
    ])

    with_hit = {h["query"] for h in hits}
    n_with = sum(1 for r in panx if r["seq_id"] in with_hit)
    pct = f"{100.0 * n_with / max(len(panx), 1):.0f}%"
    med = float(np.median([h["pident"] for h in hits])) if hits else float("nan")

    if n_with >= max(1, int(0.5 * len(panx))):
        verdict = (
            f"Most pannexin proteins ({n_with}/{len(panx)}) hit at least one innexin at e≤1e-3 "
            f"and ≥{MIN_ALN} aa (median %id {med:.1f}). That is the opposite pattern from "
            f"innexin×connexin on the parent site (no cross-family hits at the same thresholds)."
        )
    elif hits:
        verdict = (
            f"Some cross-family signal: {n_with}/{len(panx)} pannexins hit innexins "
            f"(median %id {med:.1f}). Relatedness is detectable but not uniform across the panel."
        )
    else:
        verdict = (
            "No pannexin→innexin hits at e≤1e-3 and ≥80 aa in this panel. "
            "That would be unexpected for true homologs — check panel composition and search settings."
        )

    closing = (
        "Literature treats pannexins as chordate/vertebrate relatives of innexins that no longer "
        "dock into gap junctions. This page only tests sequence detectability against the innexin "
        "panel; it does not re-derive the chordate bottleneck story from these trees alone."
    )

    plot_length(panx, inx, FIGS / "01_lengths.png")
    plot_embedding(panx, inx, FIGS / "02_kmer_space.png")
    plot_hit_identity(hits, FIGS / "03_cross_identity.png")
    plot_coverage(panx, hits, FIGS / "04_hit_coverage.png")

    stats = {
        "n_panx": str(len(panx)),
        "n_inx": str(len(inx)),
        "n_hits": str(len(hits)),
        "pct_with_hit": pct,
        "verdict": verdict,
        "closing": closing,
    }
    write_csv(OUT / "summary.csv", [stats])
    write_html(stats)
    print(f"Done → {OUT / 'index.html'}")
    print(f"  panx={len(panx)} inx={len(inx)} hits={len(hits)} with_hit={n_with} ({pct})")


if __name__ == "__main__":
    build()
