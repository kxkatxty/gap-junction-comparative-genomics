#!/usr/bin/env python3
"""Curated pannexin insights gallery (inx/cnx-style), no SynVoy.

Copies key figures from existing result pages and frames literature vs panel.
Outputs → project/results/pannexin_insights/
"""

from __future__ import annotations

import shutil
from pathlib import Path

from pipeline.common import RESULTS_DIR, read_csv_rows

OUT = RESULTS_DIR / "pannexin_insights"
FIGS = OUT / "figures"

COPIES: list[tuple[Path, str]] = [
    (RESULTS_DIR / "pannexin_path" / "figures" / "04_homology_map.png", "01_homology_map.png"),
    (RESULTS_DIR / "pannexin_path" / "figures" / "02_glycosylation_blockade.png", "02_glycosylation.png"),
    (RESULTS_DIR / "panx_vs_inx" / "figures" / "04_hit_coverage.png", "03_panx_vs_inx_coverage.png"),
    (RESULTS_DIR / "panx_vs_inx" / "figures" / "02_kmer_space.png", "04_panx_inx_kmer.png"),
    (RESULTS_DIR / "pannexin_clade_comparison" / "figures" / "02_clade_type_heatmap.png", "05_clade_type.png"),
    (RESULTS_DIR / "pannexin_clade_comparison" / "figures" / "07_bottleneck_focus.png", "06_bottleneck_copies.png"),
    (RESULTS_DIR / "pannexin_phylogeny" / "figures" / "02_iqtree_tree.png", "07_iqtree_tree.png"),
    (RESULTS_DIR / "pannexin_exon_notes" / "figures" / "01_exon_counts.png", "08_exon_counts.png"),
]


def _stats() -> dict[str, str]:
    meta = read_csv_rows(RESULTS_DIR / "pannexin_panel" / "sequence_metadata.csv")
    summary = read_csv_rows(RESULTS_DIR / "panx_vs_inx" / "summary.csv")
    copy = read_csv_rows(RESULTS_DIR / "pannexin_clade_comparison" / "species_copy_number.csv")
    n = len(meta)
    n_sp = len({r.get("organism") for r in meta})
    rescued = sum(1 for r in meta if r.get("panx_type_source") == "mmseqs_rescue")
    full = sum(1 for r in copy if r.get("has_full_triplet") == "yes")
    pct = summary[0].get("pct_with_hit", "?") if summary else "?"
    return {
        "n": str(n),
        "n_sp": str(n_sp),
        "rescued": str(rescued),
        "full": str(full),
        "pct": pct,
    }


def build() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    present: list[tuple[str, str]] = []
    for src, dest in COPIES:
        if not src.exists():
            continue
        shutil.copy2(src, FIGS / dest)
        present.append((dest, src.name))

    s = _stats()
    cards = []
    captions = {
        "01_homology_map.png": "Literature map: innexin ↔ pannexin homologous; connexins separate.",
        "02_glycosylation.png": "N-X-S/T / NGS reading: docking blocked on the pannexin branch.",
        "03_panx_vs_inx_coverage.png": f"This panel: ~{s['pct']} of pannexins hit ≥1 innexin (same thresholds as inx×cnx = 0).",
        "04_panx_inx_kmer.png": "3-mer space overlap with innexins — expected for related sequences.",
        "05_clade_type.png": "Clade × PANX1/2/3 occupancy in the UniProt panel.",
        "06_bottleneck_copies.png": "Panel copy numbers for early chordates vs sample mammals (not a genome census).",
        "07_iqtree_tree.png": "Rendered IQ-TREE on the curated panel (+ Ciona innexin outgroups).",
        "08_exon_counts.png": "Ensembl exon counts for model-species PANX genes (human / mouse / zebrafish).",
    }
    for dest, _ in present:
        cards.append(
            f"""
<article class="card">
  <img src="figures/{dest}" alt="" loading="lazy" onclick="openModal(this.src)">
  <p>{captions.get(dest, dest)}</p>
</article>"""
        )

    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pannexin insights</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600;700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root {{ --ink:#1c1917; --muted:#57534e; --ember:#c2410c; --deep:#7c2d12; --bg:#f7f3ef; }}
* {{ box-sizing:border-box }}
body {{ margin:0; font-family:Sora,sans-serif; color:var(--ink); line-height:1.55;
  background:radial-gradient(900px 420px at 0% 0%, rgba(194,65,12,.10), transparent 55%),
             linear-gradient(180deg,#faf6f1,var(--bg)); }}
.hero {{ min-height:70vh; display:grid; align-items:end; padding:clamp(1.4rem,4vw,3rem); color:#fff7ed;
  background:linear-gradient(120deg,rgba(124,45,18,.94),rgba(194,65,12,.7) 55%,rgba(15,118,110,.4)); }}
.brand {{ font-family:Fraunces,serif; font-size:clamp(2.2rem,6vw,4rem); margin:0 0 .6rem; letter-spacing:-.03em }}
.hero p {{ max-width:52ch; opacity:.93 }}
.wrap {{ width:min(1100px,calc(100% - 2rem)); margin:0 auto; padding:2rem 0 4rem }}
.stats {{ display:grid; grid-template-columns:repeat(4,1fr); gap:.75rem; margin-top:-1.6rem }}
.stat {{ background:#fff; border-radius:14px; padding:1rem; box-shadow:0 8px 22px rgba(28,25,23,.06) }}
.stat b {{ display:block; font-family:Fraunces,serif; font-size:1.6rem; color:var(--deep) }}
.stat span {{ color:var(--muted); font-size:.82rem }}
h2 {{ font-family:Fraunces,serif; font-size:1.5rem; margin:2rem 0 .6rem }}
.frame {{ background:#fff; border-left:4px solid var(--ember); border-radius:0 12px 12px 0; padding:1rem 1.1rem; color:var(--muted) }}
.grid {{ display:grid; grid-template-columns:repeat(auto-fill,minmax(300px,1fr)); gap:1rem; margin-top:1rem }}
.card {{ background:#fff; border-radius:14px; padding:.8rem; box-shadow:0 6px 18px rgba(28,25,23,.05) }}
.card img {{ width:100%; border-radius:10px; cursor:zoom-in; background:#f8fafc }}
.card p {{ margin:.55rem 0 0; font-size:.88rem; color:var(--muted) }}
.links a {{ margin-right:1rem; color:var(--ember) }}
.modal {{ display:none; position:fixed; inset:0; background:rgba(0,0,0,.88); z-index:20; align-items:center; justify-content:center; padding:1rem; cursor:zoom-out }}
.modal.open {{ display:flex }}
.modal img {{ max-width:96vw; max-height:92vh; border-radius:8px }}
@media(max-width:700px){{ .stats {{ grid-template-columns:1fr 1fr }} }}
</style>
</head>
<body>
<header class="hero">
  <div>
    <p class="brand">Pannexin insights</p>
    <p>Curated figures for the third family: homologous to innexins, not to connexins.
    Numbers and trees below come from the UniProt panel built for this thesis.</p>
  </div>
</header>
<main class="wrap">
  <div class="stats">
    <div class="stat"><b>{s['n']}</b><span>panel proteins</span></div>
    <div class="stat"><b>{s['n_sp']}</b><span>species</span></div>
    <div class="stat"><b>{s['pct']}</b><span>with ≥1 innexin hit</span></div>
    <div class="stat"><b>{s['full']}</b><span>full PANX1–3 triplets</span></div>
  </div>

  <h2>Literature vs this panel</h2>
  <div class="frame">
    <p><strong>Literature:</strong> early chordate innexin bottleneck → glycosylated remnant (pannexins) → later connexin diversification (Welzel &amp; Schuster 2022).</p>
    <p><strong>This panel:</strong> UniProt sequences, clade × type counts, panx×inx MMseqs (~{s['pct']} of pannexins hit ≥1 innexin at the same thresholds where innexin×connexin finds none),
    and an IQ-TREE on the curated set. A few generic UniProt symbols were typed by best hit to named paralogs
    ({s['rescued']} reassigned). Lancelet / lamprey copy numbers reflect what is in UniProt — not a full genome census.</p>
  </div>

  <h2>Figures</h2>
  <div class="grid">{"".join(cards)}</div>

  <p class="links" style="margin-top:2rem">
    <a href="../pannexin_path/index.html">Path</a>
    <a href="../pannexin_phylogeny/index.html">Phylogeny</a>
    <a href="../panx_vs_inx/index.html">× innexin</a>
    <a href="../pannexin_clade_comparison/index.html">Clade × type</a>
    <a href="../pannexin_exon_notes/index.html">Exon notes</a>
    <a href="../../../../project/results/phylogenetic_story/index.html">Main story</a>
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
</body>
</html>
"""
    (OUT / "index.html").write_text(page, encoding="utf-8")
    print(f"Insights → {OUT / 'index.html'} ({len(present)} figures)")


if __name__ == "__main__":
    build()
