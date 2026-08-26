#!/usr/bin/env python3
"""Annotation-gap notes for early-chordate / cartilaginous coverage in the UniProt panel.

Outputs → project/results/pannexin_curator/
"""

from __future__ import annotations

import csv
import html
from pathlib import Path

from pipeline.common import RESULTS_DIR, write_csv

OUT = RESULTS_DIR / "pannexin_curator"
STATUS_CSV = OUT / "annotation_gap_status.csv"


def write_status_table() -> None:
    rows = [
        {
            "species": "Callorhinchus milii",
            "clade": "Cartilaginous fish",
            "uniprot_pannexin": "yes (TrEMBL)",
            "action": "in curated panel (covers the skate UniProt gap)",
            "status": "in_panel",
        },
        {
            "species": "Leucoraja erinaceus",
            "clade": "Cartilaginous fish",
            "uniprot_pannexin": "none",
            "action": "no usable UniProt pannexins; elephant shark used as proxy",
            "status": "covered_by_proxy",
        },
        {
            "species": "Ciona intestinalis",
            "clade": "Tunicate / lancelet",
            "uniprot_pannexin": "none",
            "action": "UniProt shows innexins, not pannexins — used as tree outgroup only",
            "status": "outgroup_inx",
        },
        {
            "species": "Branchiostoma floridae",
            "clade": "Tunicate / lancelet",
            "uniprot_pannexin": "LOC / pannexin-like labels",
            "action": "in curated panel (two proteins); orthology still tentative",
            "status": "in_panel_weak",
        },
        {
            "species": "Petromyzon marinus",
            "clade": "Jawless vertebrate",
            "uniprot_pannexin": "LOC labels",
            "action": "in curated panel (one protein); orthology still tentative",
            "status": "in_panel_weak",
        },
    ]
    write_csv(STATUS_CSV, rows)


def write_html() -> None:
    status_rows = list(csv.DictReader(STATUS_CSV.open(encoding="utf-8"))) if STATUS_CSV.exists() else []
    status_html = "".join(
        f"<tr><td>{html.escape(r['species'])}</td><td>{html.escape(r['clade'])}</td>"
        f"<td>{html.escape(r['uniprot_pannexin'])}</td><td>{html.escape(r['status'])}</td>"
        f"<td>{html.escape(r['action'])}</td></tr>"
        for r in status_rows
    )
    page = f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pannexin annotation gaps</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root {{ --ink:#12202a; --muted:#4d6570; --ember:#c2410c; --line:rgba(18,32,42,.12); }}
body {{ margin:0; font-family:Sora,sans-serif; color:var(--ink); background:#f7f3ef; line-height:1.55 }}
.wrap {{ width:min(920px, calc(100% - 2rem)); margin:0 auto; padding:2rem 0 4rem }}
h1 {{ font-family:Fraunces,serif; margin:0 0 .4rem }}
.lede {{ color:var(--muted); max-width:60ch }}
.take {{ border-left:3px solid var(--ember); background:rgba(194,65,12,.07); padding:.75rem .9rem; border-radius:0 12px 12px 0 }}
.table-wrap {{ overflow:auto; background:#fff; border:1px solid var(--line); border-radius:12px; margin:1rem 0 }}
table {{ border-collapse:collapse; width:100%; font-size:.82rem }}
th,td {{ padding:.45rem .55rem; border-bottom:1px solid var(--line); text-align:left }}
th {{ background:rgba(194,65,12,.08) }}
.footer a {{ margin-right:.85rem; color:var(--ember) }}
</style></head><body>
<main class="wrap">
  <h1>Pannexin annotation gaps</h1>
  <p class="lede">Where UniProt coverage is thin for the chordate story — and how the curated panel
  handles it. This page is a status table, not a genome search.</p>
  <p class="take"><strong>In practice.</strong> The cartilaginous gap is covered with
  <em>Callorhinchus milii</em>. Tunicates contribute innexins (tree outgroups), not pannexins.
  Lancelet and lamprey entries in the panel still carry weak LOC labels.</p>
  <div class="table-wrap"><table>
    <thead><tr><th>species</th><th>clade</th><th>UniProt panx</th><th>status</th><th>action</th></tr></thead>
    <tbody>{status_html}</tbody>
  </table></div>
  <p class="footer">
    <a href="../pannexin_annotation_status/index.html">Annotation status</a>
    <a href="../pannexin_clade_comparison/index.html">Clade × type</a>
    <a href="../panx_vs_inx/index.html">Pannexin × innexin</a>
    <a href="../pannexin_phylogeny/index.html">Phylogeny</a>
    <a href="../pannexin_path/index.html">Path</a>
  </p>
</main></body></html>
"""
    (OUT / "index.html").write_text(page, encoding="utf-8")


def build() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    write_status_table()
    write_html()
    print(f"Done → {OUT / 'index.html'}")


if __name__ == "__main__":
    build()
