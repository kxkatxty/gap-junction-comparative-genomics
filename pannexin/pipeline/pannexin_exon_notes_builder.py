#!/usr/bin/env python3
"""Light exon notes for model-species pannexins via Ensembl REST (no SynVoy).

Fetches coding-exon counts for PANX1/2/3 in human, mouse, zebrafish when network
is available. Documents that this is not a genome-wide SynVoy architecture survey.

Outputs → project/results/pannexin_exon_notes/
"""

from __future__ import annotations

import html
import json
import urllib.error
import urllib.request
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pipeline.common import RESULTS_DIR, TYPE_COLORS, write_csv

OUT = RESULTS_DIR / "pannexin_exon_notes"
FIGS = OUT / "figures"
CACHE = OUT / "ensembl_exon_counts.csv"

# Ensembl gene symbols / species
QUERIES = [
    ("Homo sapiens", "human", "PANX1"),
    ("Homo sapiens", "human", "PANX2"),
    ("Homo sapiens", "human", "PANX3"),
    ("Mus musculus", "mouse", "Panx1"),
    ("Mus musculus", "mouse", "Panx2"),
    ("Mus musculus", "mouse", "Panx3"),
    ("Danio rerio", "zebrafish", "panx1a"),
    ("Danio rerio", "zebrafish", "panx1b"),
    ("Danio rerio", "zebrafish", "panx2"),
    ("Danio rerio", "zebrafish", "panx3"),
]

UA = "BachelorCursor-pannexin-exons/1.0"


def _get_json(url: str) -> dict | list | None:
    req = urllib.request.Request(url, headers={"User-Agent": UA, "Content-Type": "application/json"})
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            return json.loads(resp.read().decode("utf-8"))
    except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError):
        return None


def fetch_exon_count(species: str, symbol: str) -> dict[str, str]:
    lookup = _get_json(f"https://rest.ensembl.org/lookup/symbol/{species}/{symbol}?content-type=application/json")
    if not isinstance(lookup, dict) or "id" not in lookup:
        return {
            "organism": "",
            "ensembl_species": species,
            "symbol": symbol,
            "gene_id": "",
            "n_coding_exons": "",
            "status": "not_found",
        }
    gene_id = lookup["id"]
    # transcript overlap exons
    overlap = _get_json(
        f"https://rest.ensembl.org/overlap/id/{gene_id}?feature=exon&content-type=application/json"
    )
    n_exons = ""
    status = "ok"
    if isinstance(overlap, list):
        # unique exon ids
        ids = {e.get("id") for e in overlap if isinstance(e, dict) and e.get("id")}
        n_exons = str(len(ids)) if ids else "0"
    else:
        status = "overlap_failed"
    return {
        "organism": "",
        "ensembl_species": species,
        "symbol": symbol,
        "gene_id": gene_id,
        "n_coding_exons": n_exons,
        "status": status,
    }


def collect_rows() -> list[dict[str, str]]:
    rows = []
    for organism, ensembl_sp, symbol in QUERIES:
        row = fetch_exon_count(ensembl_sp, symbol)
        row["organism"] = organism
        row["panx_type"] = (
            "PANX1"
            if "1" in symbol.upper()
            else "PANX2"
            if "2" in symbol.upper()
            else "PANX3"
            if "3" in symbol.upper()
            else "other/unknown"
        )
        rows.append(row)
        print(f"  {organism} {symbol}: exons={row['n_coding_exons'] or '?'} ({row['status']})")
    return rows


def plot_exons(rows: list[dict[str, str]], path: Path) -> None:
    usable = [r for r in rows if r.get("n_coding_exons") and r["n_coding_exons"].isdigit()]
    if not usable:
        return
    labels = [f"{r['organism'].split()[0][0]}. {r['symbol']}" for r in usable]
    vals = [int(r["n_coding_exons"]) for r in usable]
    colors = [TYPE_COLORS.get(r["panx_type"], "#94A3B8") for r in usable]
    fig, ax = plt.subplots(figsize=(9.5, 4.8))
    x = np.arange(len(vals))
    ax.bar(x, vals, color=colors, width=0.72)
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=30, ha="right", fontsize=8)
    ax.set_ylabel("Ensembl exon features (gene overlap)")
    ax.set_title("Model-species PANX loci — light exon counts (Ensembl REST)")
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, dpi=180, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def write_html(rows: list[dict[str, str]]) -> None:
    table = "".join(
        "<tr>"
        f"<td>{html.escape(r['organism'])}</td>"
        f"<td>{html.escape(r['symbol'])}</td>"
        f"<td>{html.escape(r['panx_type'])}</td>"
        f"<td>{html.escape(r['n_coding_exons'] or '—')}</td>"
        f"<td>{html.escape(r['gene_id'])}</td>"
        f"<td>{html.escape(r['status'])}</td>"
        "</tr>"
        for r in rows
    )
    fig_block = ""
    if (FIGS / "01_exon_counts.png").exists():
        fig_block = """
  <div class="fig"><img src="figures/01_exon_counts.png" alt="exon counts">
    <p class="cap">Ensembl exon-feature counts for selected PANX genes.</p></div>"""
    page = f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pannexin exon notes</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root {{ --ink:#12202a; --muted:#4d6570; --ember:#c2410c; --line:rgba(18,32,42,.12); }}
body {{ margin:0; font-family:Sora,sans-serif; color:var(--ink); background:#f7f3ef; line-height:1.55 }}
.wrap {{ width:min(920px, calc(100% - 2rem)); margin:0 auto; padding:2rem 0 4rem }}
h1 {{ font-family:Fraunces,serif; margin:0 0 .4rem }}
.lede {{ color:var(--muted); max-width:64ch }}
.take {{ border-left:3px solid var(--ember); background:rgba(194,65,12,.07); padding:.75rem .9rem; border-radius:0 12px 12px 0 }}
.fig {{ background:#fff; border:1px solid var(--line); border-radius:14px; padding:.85rem; margin:1rem 0 }}
.fig img {{ width:100%; display:block; border-radius:8px }}
.cap {{ font-size:.88rem; color:var(--muted); margin:.5rem 0 0 }}
.table-wrap {{ overflow:auto; background:#fff; border:1px solid var(--line); border-radius:12px; margin:1rem 0 }}
table {{ border-collapse:collapse; width:100%; font-size:.82rem }}
th,td {{ padding:.45rem .55rem; border-bottom:1px solid var(--line); text-align:left }}
th {{ background:rgba(194,65,12,.08) }}
.footer a {{ margin-right:.85rem; color:var(--ember) }}
</style></head><body>
<main class="wrap">
  <h1>Pannexin exon notes</h1>
  <p class="lede">Ensembl REST lookup for model-species PANX genes — a compact counterpart to the
  exon-architecture panels on the innexin/connexin sites, limited to human, mouse and zebrafish.</p>
  <p class="take"><strong>Scope.</strong> Counts are Ensembl exon features overlapping each gene.
  They support the simple claim that mammalian (and zebrafish) PANX loci are multi-exon;
  they are not a curated coding-exon gold standard across the full panel.</p>
  {fig_block}
  <div class="table-wrap"><table>
    <thead><tr><th>organism</th><th>symbol</th><th>type</th><th>exons</th><th>Ensembl id</th><th>status</th></tr></thead>
    <tbody>{table}</tbody>
  </table></div>
  <p class="footer">
    <a href="../pannexin_insights/index.html">Insights</a>
    <a href="../pannexin_path/index.html">Path</a>
    <a href="../pannexin_clade_comparison/index.html">Clade × type</a>
  </p>
</main></body></html>
"""
    (OUT / "index.html").write_text(page, encoding="utf-8")


def build(*, reuse: bool = False) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    if reuse and CACHE.exists():
        from pipeline.common import read_csv_rows

        rows = read_csv_rows(CACHE)
        print(f"Reusing {CACHE}")
    else:
        print("Fetching Ensembl exon counts (model species)…")
        rows = collect_rows()
        write_csv(CACHE, rows)
    plot_exons(rows, FIGS / "01_exon_counts.png")
    write_html(rows)
    print(f"Done → {OUT / 'index.html'}")


if __name__ == "__main__":
    import sys

    build(reuse="--reuse" in sys.argv)
