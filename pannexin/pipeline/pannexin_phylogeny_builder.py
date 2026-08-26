#!/usr/bin/env python3
"""Simple PANX1/2/3 phylogeny (MAFFT + IQ-TREE) for the curated panel.

Optionally adds Ciona innexin outgroups. No connexins.

Outputs → project/results/pannexin_phylogeny/
"""

from __future__ import annotations

import html
import re
import shutil
import subprocess
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, Patch

from pipeline.common import PROJECT_ROOT, RESULTS_DIR, TYPE_COLORS, read_csv_rows
from pipeline.limits import n_threads, run_heavy_phylogeny
from pipeline.newick_plot import iter_edges, iter_tips, ladderize, layout, parse_newick
from pipeline.pannexin_panel_builder import FASTA as PANEL_FASTA
from pipeline.pannexin_panel_builder import META as PANEL_META
from pipeline.pannexin_panel_builder import build as build_panel

OUT = RESULTS_DIR / "pannexin_phylogeny"
FIGS = OUT / "figures"
WORK = OUT / "work"
ALIGNED = OUT / "pannexin_panel.aln.fasta"
TREE = OUT / "pannexin_panel.treefile"
OUTGROUP_DIR = PROJECT_ROOT / "project" / "data" / "references" / "outgroups" / "ciona_innexins"

SYNVOY_BIN = Path.home() / "miniconda3" / "envs" / "synvoy_env" / "bin"
MAFFT = str(SYNVOY_BIN / "mafft") if (SYNVOY_BIN / "mafft").exists() else (shutil.which("mafft") or "mafft")
IQTREE = str(SYNVOY_BIN / "iqtree") if (SYNVOY_BIN / "iqtree").exists() else (
    shutil.which("iqtree2") or shutil.which("iqtree") or "iqtree"
)


def _clean_id(seq_id: str) -> str:
    # IQ-TREE / Newick-safe-ish labels
    return re.sub(r"[^A-Za-z0-9_.|-]", "_", seq_id)[:80]


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    if not path.exists():
        return out
    hdr, chunks = None, []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(chunks)))
            hdr = line[1:].strip().split()[0]
            chunks = []
        else:
            chunks.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(chunks)))
    return out


def build_input_fasta(dest: Path) -> tuple[int, int]:
    build_panel()
    rows = _read_fasta(PANEL_FASTA)
    meta = {r["seq_id"]: r for r in read_csv_rows(PANEL_META)}
    # Prefer typed paralogs for a readable tree; keep a few other/unknown chordates
    typed = []
    other = []
    for sid, seq in rows:
        m = meta.get(sid, {})
        if m.get("panx_type") in {"PANX1", "PANX2", "PANX3"}:
            typed.append((sid, seq, m.get("panx_type", "?")))
        else:
            other.append((sid, seq, m.get("panx_type", "other/unknown")))
    # Cap other/unknown to avoid clutter (lamprey/lancelet)
    other = other[:6]

    outgroup = []
    for path in sorted(OUTGROUP_DIR.glob("*.fasta")) if OUTGROUP_DIR.exists() else []:
        for hdr, seq in _read_fasta(path):
            if len(seq) >= 200:
                oid = _clean_id(f"out|Ciona|{path.stem}")
                outgroup.append((oid, seq, "outgroup_inx"))

    dest.parent.mkdir(parents=True, exist_ok=True)
    with dest.open("w", encoding="utf-8") as handle:
        for sid, seq, _ in typed + other + outgroup:
            handle.write(f">{_clean_id(sid)}\n")
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")
    return len(typed) + len(other), len(outgroup)


def run_mafft(inp: Path, out: Path) -> None:
    log = out.with_suffix(".mafft.log")
    threads = str(n_threads(1))
    # FFT-NS-2 is much lighter than L-INS-i from --auto on ~50 proteins.
    cmd = [MAFFT, "--retree", "2", "--maxiterate", "0", "--thread", threads, str(inp)]
    with log.open("w", encoding="utf-8") as err:
        proc = subprocess.run(cmd, check=True, stdout=subprocess.PIPE, stderr=err, text=True)
    out.write_text(proc.stdout, encoding="utf-8")
    if out.stat().st_size < 100:
        raise RuntimeError(f"MAFFT produced empty alignment; see {log}")


def run_iqtree(aln: Path, prefix: Path) -> Path:
    prefix.parent.mkdir(parents=True, exist_ok=True)
    for ext in (".treefile", ".iqtree", ".log", ".ckp.gz", ".bionj", ".mldist", ".contree", ".splits.nex"):
        p = Path(str(prefix) + ext)
        if p.exists():
            p.unlink()
    log = Path(str(prefix) + ".run.log")
    threads = str(n_threads(1))
    # Light default: fixed model, no UFBoot. Heavy mode adds -bb 200 only.
    cmd = [
        IQTREE,
        "-s",
        str(aln),
        "-m",
        "LG+G4",
        "-nt",
        threads,
        "-pre",
        str(prefix),
    ]
    if run_heavy_phylogeny():
        cmd[cmd.index("-m") + 1] = "LG+I+G4"
        cmd.extend(["-bb", "200"])
    with log.open("w", encoding="utf-8") as err:
        subprocess.run(cmd, check=True, stdout=err, stderr=subprocess.STDOUT)
    treefile = Path(str(prefix) + ".treefile")
    if not treefile.exists():
        raise FileNotFoundError(treefile)
    shutil.copy2(treefile, TREE)
    return TREE


def _tip_type(name: str | None, meta: dict[str, dict]) -> str:
    n = name or ""
    if n.upper().startswith("OUT|") or n.upper().startswith("OUTGROUP"):
        return "outgroup_inx"
    if n in meta:
        return meta[n].get("panx_type") or "other/unknown"
    up = n.upper()
    for key in ("PANX1", "PANX2", "PANX3"):
        if key in up.replace("-", "").replace("_", ""):
            return key
    return "other/unknown"


def _short_label(name: str | None, meta: dict[str, dict]) -> str:
    n = name or ""
    if n.upper().startswith("OUT|"):
        parts = n.split("|")
        return f"out:{parts[-1][:16]}" if len(parts) > 1 else "outgroup"
    m = meta.get(n, {})
    org = (m.get("organism") or "").replace(" ", "_")
    if not org and "|" in n:
        bits = n.split("|")
        org = bits[1] if len(bits) > 1 else ""
        gene = bits[2] if len(bits) > 2 else ""
    else:
        gene = m.get("gene") or ""
    org_short = org.split("_")[0][:10] if org else "?"
    gene_short = (gene or "?")[:12]
    typ = _tip_type(n, meta)
    prefix = {"PANX1": "1", "PANX2": "2", "PANX3": "3"}.get(typ, "?")
    return f"{prefix}:{org_short}/{gene_short}"


def plot_iqtree_tree(meta_rows: list[dict], path: Path) -> None:
    """Render the actual IQ-TREE Newick as a rectangular phylogram."""
    if not TREE.exists():
        raise FileNotFoundError(TREE)
    meta = {r["seq_id"]: r for r in meta_rows}
    root = parse_newick(TREE.read_text(encoding="utf-8"))
    ladderize(root)
    layout(root)
    tips = list(iter_tips(root))
    n_tips = len(tips)
    fig_h = max(8.0, 0.28 * n_tips + 1.5)
    fig, ax = plt.subplots(figsize=(11.5, fig_h))
    colors = {
        "PANX1": TYPE_COLORS["PANX1"],
        "PANX2": TYPE_COLORS["PANX2"],
        "PANX3": TYPE_COLORS["PANX3"],
        "other/unknown": TYPE_COLORS["other/unknown"],
        "outgroup_inx": "#64748B",
    }
    for parent, child in iter_edges(root):
        # elbow: horizontal then vertical
        ax.plot([parent.x, parent.x], [parent.y, child.y], color="#94A3B8", lw=0.9, solid_capstyle="round")
        ax.plot([parent.x, child.x], [child.y, child.y], color="#94A3B8", lw=0.9, solid_capstyle="round")
    xmax = max((t.x for t in tips), default=1.0)
    for tip in tips:
        typ = _tip_type(tip.name, meta)
        col = colors.get(typ, "#94A3B8")
        ax.plot([tip.x], [tip.y], "o", color=col, markersize=4.5, zorder=3)
        ax.text(
            tip.x + xmax * 0.012,
            tip.y,
            _short_label(tip.name, meta),
            va="center",
            ha="left",
            fontsize=7.2,
            color=col,
            fontweight="600",
        )
    ax.set_ylim(-1, n_tips)
    ax.set_xlim(-xmax * 0.02, xmax * 1.38)
    ax.invert_yaxis()
    ax.set_yticks([])
    ax.set_xlabel("Substitutions / site (IQ-TREE branch length)")
    ax.set_title("Pannexin panel phylogeny (IQ-TREE) — tips colored by PANX type", loc="left", fontsize=12)
    handles = [
        Patch(color=colors["PANX1"], label="PANX1"),
        Patch(color=colors["PANX2"], label="PANX2"),
        Patch(color=colors["PANX3"], label="PANX3"),
        Patch(color=colors["other/unknown"], label="other / weak"),
        Patch(color=colors["outgroup_inx"], label="Ciona innexin outgroup"),
    ]
    ax.legend(handles=handles, loc="lower left", frameon=False, fontsize=8)
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["left"].set_visible(False)
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def plot_schematic(meta_rows: list[dict], n_out: int, path: Path) -> None:
    """Small expectation card kept as a companion figure."""
    counts = {"PANX1": 0, "PANX2": 0, "PANX3": 0, "other/unknown": 0}
    for r in meta_rows:
        t = r.get("panx_type", "other/unknown")
        if t in counts:
            counts[t] += 1
        else:
            counts["other/unknown"] += 1
    fig, ax = plt.subplots(figsize=(8.2, 3.6))
    ax.set_xlim(0, 10)
    ax.set_ylim(0, 5)
    ax.axis("off")
    ax.set_title("Reading aid: expected PANX1 / PANX2 / PANX3 blocks after WGD", loc="left", fontsize=11)
    ax.plot([1.2, 2.2], [2.5, 2.5], color="#334155", lw=2)
    ypos = {"PANX1": 3.8, "PANX2": 2.5, "PANX3": 1.2}
    for name, y in ypos.items():
        ax.plot([2.2, 3.4], [2.5, y], color="#334155", lw=1.5)
        ax.add_patch(
            FancyBboxPatch(
                (3.5, y - 0.4),
                5.8,
                0.8,
                boxstyle="round,pad=0.02,rounding_size=0.2",
                facecolor=TYPE_COLORS.get(name, "#94A3B8"),
                edgecolor="#1e293b",
                alpha=0.85,
            )
        )
        ax.text(
            6.4,
            y,
            f"{name}  (n={counts.get(name, 0)})",
            ha="center",
            va="center",
            color="white",
            fontweight="bold",
            fontsize=11,
        )
    ax.text(
        5,
        0.35,
        f"Compare to the rendered IQ-TREE figure · Ciona innexin outgroups: {n_out}",
        ha="center",
        fontsize=8.5,
        color="#64748b",
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_html(n_seq: int, n_out: int, model_line: str) -> None:
    page = f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pannexin phylogeny</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root {{ --ink:#12202a; --muted:#4d6570; --ember:#c2410c; --line:rgba(18,32,42,.12); }}
body {{ margin:0; font-family:Sora,sans-serif; color:var(--ink); background:#f7f3ef;
  background-image:radial-gradient(900px 420px at 90% 0%, rgba(194,65,12,.10), transparent 55%); }}
.wrap {{ width:min(980px, calc(100% - 2rem)); margin:0 auto; padding:2rem 0 4rem }}
h1 {{ font-family:Fraunces,serif; margin:0 0 .4rem }}
.lede {{ color:var(--muted); max-width:68ch }}
.fig {{ background:#fff; border:1px solid var(--line); border-radius:14px; padding:.85rem; margin:1rem 0 }}
.fig img {{ width:100%; display:block; border-radius:8px; cursor:zoom-in }}
.cap {{ font-size:.88rem; color:var(--muted); margin:.5rem 0 0 }}
.footer a {{ margin-right:.85rem; color:var(--ember) }}
.take {{ border-left:3px solid var(--ember); background:rgba(194,65,12,.07); padding:.75rem .9rem; border-radius:0 12px 12px 0 }}
.note {{ border-left:3px solid #0e7490; background:rgba(14,116,144,.07); padding:.75rem .9rem; border-radius:0 12px 12px 0; margin:1rem 0; color:var(--muted); font-size:.92rem }}
.modal {{ display:none; position:fixed; inset:0; background:rgba(0,0,0,.88); z-index:20; align-items:center; justify-content:center; padding:1rem; cursor:zoom-out }}
.modal.open {{ display:flex }}
.modal img {{ max-width:96vw; max-height:92vh; border-radius:8px }}
</style></head><body>
<main class="wrap">
  <h1>Pannexin phylogeny</h1>
  <p class="lede">MAFFT alignment + IQ-TREE on the curated vertebrate/chordate
  pannexin panel (<strong>{n_seq}</strong> sequences + <strong>{n_out}</strong> Ciona innexin outgroup sequences).
  Connexins are not included. The figure below is the <strong>rendered treefile</strong>, not a cartoon.</p>
  <p class="take"><strong>How to read it.</strong> Tip colours follow panel types (PANX1 / PANX2 / PANX3;
  a few generic UniProt symbols were typed by best MMseqs hit to named paralogs).
  Expect separate vertebrate paralog blocks after whole-genome duplications; teleosts can show
  extra PANX1 copies. Early-chordate tips are sparse because few curated UniProt entries exist —
  that is annotation occupancy, not a claim about every gene in every genome.
  Substitution model: {html.escape(model_line or "LG+G4")}.</p>
  <div class="fig"><img src="figures/02_iqtree_tree.png" alt="IQ-TREE pannexin phylogeny" onclick="openModal(this.src)">
    <p class="cap">IQ-TREE phylogram from the curated alignment. Labels are shortened (type:Genus/gene).</p></div>
  <div class="fig"><img src="figures/01_expected_clades.png" alt="expected clades reading aid">
    <p class="cap">Reading aid: expected PANX1–3 blocks after vertebrate genome duplications.</p></div>
  <p class="note"><strong>Literature vs this panel.</strong> The early-chordate bottleneck → pannexin → connexin
  sequence is the published evolutionary framing (Welzel &amp; Schuster). This tree asks a narrower question:
  do the curated sequences fall into PANX1 / PANX2 / PANX3-like groups?</p>
  <p class="footer">
    <a href="pannexin_panel.treefile">Download treefile</a>
    <a href="pannexin_panel.aln.fasta">Alignment</a>
    <a href="../pannexin_insights/index.html">Insights</a>
    <a href="../panx_vs_inx/index.html">Pannexin × innexin</a>
    <a href="../pannexin_path/index.html">Path</a>
  </p>
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
</body></html>
"""
    (OUT / "index.html").write_text(page, encoding="utf-8")


def _n_outgroups_in_tree() -> int:
    if not TREE.exists():
        return 0
    text = TREE.read_text(encoding="utf-8")
    return len(re.findall(r"out\|", text, flags=re.I))


def build() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    WORK.mkdir(parents=True, exist_ok=True)

    # Reuse an existing treefile when present — avoids re-running IQ-TREE on every path rebuild.
    if TREE.exists() and ALIGNED.exists() and not run_heavy_phylogeny():
        meta = read_csv_rows(PANEL_META)
        n_out = _n_outgroups_in_tree()
        plot_iqtree_tree(meta, FIGS / "02_iqtree_tree.png")
        plot_schematic(meta, n_out, FIGS / "01_expected_clades.png")
        write_html(len(meta), n_out, "existing alignment and tree (no recompute)")
        print(f"Reused existing tree → {TREE}")
        print(f"Done → {OUT / 'index.html'}")
        return

    raw = WORK / "input.fasta"
    n_seq, n_out = build_input_fasta(raw)
    print(f"Phylogeny input: {n_seq} pannexins + {n_out} outgroups (threads={n_threads(1)})")
    run_mafft(raw, ALIGNED)
    prefix = WORK / "pannexin_panel"
    run_iqtree(ALIGNED, prefix)
    report = Path(str(prefix) + ".iqtree")
    model_line = "LG+G4 (light)" if not run_heavy_phylogeny() else "LG+I+G4 + UFBoot 200"
    if report.exists():
        for line in report.read_text(encoding="utf-8", errors="replace").splitlines():
            if "Best-fit model" in line or "Model of substitution" in line:
                model_line = line.strip()
                break
        shutil.copy2(report, OUT / "pannexin_panel.iqtree")
    meta = read_csv_rows(PANEL_META)
    plot_iqtree_tree(meta, FIGS / "02_iqtree_tree.png")
    plot_schematic(meta, n_out, FIGS / "01_expected_clades.png")
    write_html(n_seq, n_out, model_line)
    print(f"Done → {OUT / 'index.html'}")
    print(f"Tree → {TREE}")


if __name__ == "__main__":
    build()
