#!/usr/bin/env python3
"""Build a showcase website for curator-recovered innexin species."""

from __future__ import annotations

import csv
import html
from collections import defaultdict
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SUMMARY = ROOT / "project/results/gene_curator_probe/innexin_search_summary.csv"
CAND = ROOT / "project/metadata/gap_junction_candidates_master.csv"
OUT_DIR = ROOT / "project/results/new_species_gallery"
OUT_HTML = OUT_DIR / "index.html"

CLADES = {
    "Abra_alba": ("Mollusca", "bivalve"),
    "Abra_segmentum": ("Mollusca", "bivalve"),
    "Acanthocardia_echinata": ("Mollusca", "bivalve"),
    "Acanthochitona_discrepans": ("Mollusca", "chiton"),
    "Acartia_tonsa": ("Arthropoda", "copepod"),
    "Adineta_ricciae": ("Rotifera", "bdelloid"),
    "Adineta_steineri": ("Rotifera", "bdelloid"),
    "Adineta_vaga": ("Rotifera", "bdelloid"),
    "Aegaeobuthus_cyprius": ("Arthropoda", "scorpion"),
    "Aelurillus_cypriotus": ("Arthropoda", "spider"),
    "Agelena_orientalis": ("Arthropoda", "spider"),
    "Amaurobius_ferox": ("Arthropoda", "spider"),
    "Brachionus_calyciflorus": ("Rotifera", "monogonont"),
    "Brachionus_koreanus": ("Rotifera", "monogonont"),
    "Brachionus_manjavacas": ("Rotifera", "monogonont"),
    "Dryodora_glandiformis": ("Ctenophora", "ctenophore"),
    "Eriophyidae_sp.": ("Arthropoda", "mite"),
    "Rotaria_macrura": ("Rotifera", "bdelloid"),
    "Streblospio_benedicti": ("Annelida", "polychaete"),
    "Triops_cancriformis": ("Arthropoda", "branchiopod"),
}

CLADES_ORDER = ["Rotifera", "Mollusca", "Arthropoda", "Annelida", "Ctenophora"]

NOTES = {
    "Rotaria_macrura": "Rescued: prior search had candidates but zero accepted loci.",
    "Brachionus_manjavacas": "Rescued: prior accepted=0; curator recovered present models.",
    "Abra_alba": "Rescued: prior accepted=0; now a clear multi-locus molluscan hit.",
    "Agelena_orientalis": "New spider genome with strong UNC-9/Inx-like identity (~56%).",
    "Amaurobius_ferox": "New spider expansion — 11 present loci.",
    "Aelurillus_cypriotus": "New jumping-spider genome; large-genome round.",
    "Acanthocardia_echinata": "Highest molluscan yield in this batch (14 present loci).",
    "Streblospio_benedicti": "First polychaete with a clear multi-locus curator recovery.",
    "Adineta_steineri": "Highest identity among new rotifers (best id ~73%).",
    "Acartia_tonsa": "Copepod with strong best-hit identity (~65%).",
}


def _is_new(row: dict) -> bool:
    prior = str(row.get("prior_discovery") or "").strip()
    if prior in ("", "not_run"):
        return True
    try:
        return int(float(row.get("prior_accepted") or 0)) == 0
    except (TypeError, ValueError):
        return False


def _display_name(species_dir: str) -> str:
    return species_dir.replace("_", " ")


def _load_best_hits() -> dict[str, str]:
    hits: dict[str, str] = {}
    if not CAND.exists():
        return hits
    with CAND.open(newline="") as fh:
        for row in csv.DictReader(fh):
            if row.get("discovery_batch") != "gene_curator_probe":
                continue
            if row.get("rank_category") != "curator_present_innexin":
                continue
            sp = row.get("species_dir") or ""
            hit = (row.get("best_hit") or "").strip()
            if not sp or not hit:
                continue
            # keep first / shortest readable token
            token = hit.split()[0]
            prev = hits.get(sp)
            if prev is None:
                hits[sp] = token
            elif token not in prev:
                hits[sp] = f"{prev}; {token}" if len(prev) < 80 else prev
    return hits


def _load_species() -> list[dict]:
    rows = list(csv.DictReader(SUMMARY.open(newline="")))
    best_hits = _load_best_hits()
    out = []
    for r in rows:
        present = int(r.get("present_loci") or 0)
        if present <= 0:
            continue
        sp = r["species"]
        clade, habit = CLADES.get(sp, ("Other", ""))
        new = _is_new(r)
        prior = str(r.get("prior_discovery") or "")
        if prior == "not_run":
            status = "new"
            status_label = "New"
        elif new:
            status = "rescued"
            status_label = "Rescued"
        else:
            status = "expanded"
            status_label = "Expanded"
        out.append(
            {
                "species": sp,
                "name": _display_name(sp),
                "clade": clade,
                "habit": habit,
                "present": present,
                "fragmentary": int(r.get("fragmentary_loci") or 0),
                "weak": int(r.get("weak_loci") or 0),
                "best_id": float(r.get("best_identity") or 0),
                "best_aa": int(float(r.get("best_aa") or 0)),
                "best_locus": r.get("best_locus") or "",
                "prior": prior,
                "status": status,
                "status_label": status_label,
                "is_new": new,
                "note": NOTES.get(sp, ""),
                "best_hit": best_hits.get(sp, ""),
            }
        )
    out.sort(key=lambda x: (0 if x["is_new"] else 1, -x["present"], x["name"]))
    return out


def _bar(pct: float) -> str:
    w = max(4, min(100, round(pct * 100)))
    return f'<span class="idbar"><i style="width:{w}%"></i></span>'


def _species_card(s: dict) -> str:
    note = (
        f'<p class="note">{html.escape(s["note"])}</p>' if s["note"] else ""
    )
    hit = (
        f'<div class="meta"><span>Best hit</span><b>{html.escape(s["best_hit"])}</b></div>'
        if s["best_hit"]
        else ""
    )
    habit = f' · {html.escape(s["habit"])}' if s["habit"] else ""
    return f"""
<article class="species {s['status']}" data-clade="{html.escape(s['clade'])}" data-status="{s['status']}" data-name="{html.escape(s['name'].lower())}">
  <header>
    <span class="pill {s['status']}">{s['status_label']}</span>
    <span class="clade-tag">{html.escape(s['clade'])}{habit}</span>
  </header>
  <h3><em>{html.escape(s['name'])}</em></h3>
  <div class="metrics">
    <div class="metric"><b>{s['present']}</b><span>present loci</span></div>
    <div class="metric"><b>{s['best_aa']}</b><span>best length (aa)</span></div>
    <div class="metric"><b>{s['best_id']:.0%}</b><span>best identity</span></div>
  </div>
  <div class="idrow"><span>Identity</span>{_bar(s['best_id'])}<b>{s['best_id']:.1%}</b></div>
  {hit}
  <div class="meta"><span>Best locus</span><code>{html.escape(s['best_locus'] or '—')}</code></div>
  <div class="meta"><span>Prior discovery</span><b>{html.escape(s['prior'] or '—')}</b></div>
  {note}
</article>"""


def build_html(species: list[dict]) -> str:
    new = [s for s in species if s["is_new"]]
    expanded = [s for s in species if not s["is_new"]]
    n_new = len(new)
    n_loci_new = sum(s["present"] for s in new)
    n_present = len(species)
    n_loci = sum(s["present"] for s in species)
    n_clades = len({s["clade"] for s in new})

    by_clade: dict[str, list[dict]] = defaultdict(list)
    for s in new:
        by_clade[s["clade"]].append(s)

    clade_sections = []
    for clade in CLADES_ORDER:
        items = by_clade.get(clade) or []
        if not items:
            continue
        cards = "\n".join(_species_card(s) for s in items)
        clade_sections.append(
            f"""
<section class="clade-block" id="clade-{html.escape(clade.lower())}">
  <div class="clade-head">
    <h2>{html.escape(clade)}</h2>
    <p>{len(items)} new/rescued species · {sum(s['present'] for s in items)} present loci</p>
  </div>
  <div class="grid">{cards}</div>
</section>"""
        )

    expanded_cards = "\n".join(_species_card(s) for s in expanded)

    return f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Curator Recovery — New Innexin Species</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link rel="preconnect" href="https://fonts.gstatic.com" crossorigin>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,500;9..144,700&family=Sora:wght@400;500;600;700&display=swap" rel="stylesheet">
<style>
:root {{
  --ink: #12232e;
  --muted: #4d6573;
  --foam: #f3f8f6;
  --sea: #0f6a6a;
  --sea-deep: #0a3d44;
  --sand: #d9a15b;
  --coral: #c45c3e;
  --mint: #2f9e7b;
  --card: rgba(255,255,255,.88);
  --line: rgba(18,35,46,.1);
}}
* {{ box-sizing: border-box; }}
html {{ scroll-behavior: smooth; }}
body {{
  margin: 0;
  color: var(--ink);
  font-family: "Sora", sans-serif;
  background:
    radial-gradient(1200px 600px at 10% -10%, rgba(47,158,123,.18), transparent 55%),
    radial-gradient(900px 500px at 90% 0%, rgba(217,161,91,.16), transparent 50%),
    linear-gradient(180deg, #e8f2ef 0%, var(--foam) 35%, #eef3f1 100%);
  line-height: 1.55;
}}
a {{ color: var(--sea); }}
.hero {{
  position: relative;
  min-height: 92vh;
  display: grid;
  align-items: end;
  padding: clamp(1.5rem, 4vw, 3.5rem);
  color: #f7fffc;
  overflow: hidden;
  background:
    linear-gradient(115deg, rgba(10,61,68,.92) 0%, rgba(15,106,106,.78) 48%, rgba(196,92,62,.55) 100%),
    url("data:image/svg+xml,%3Csvg xmlns='http://www.w3.org/2000/svg' width='160' height='160' viewBox='0 0 160 160'%3E%3Cg fill='none' stroke='%23ffffff' stroke-opacity='0.08' stroke-width='1'%3E%3Cpath d='M0 40h160M0 80h160M0 120h160M40 0v160M80 0v160M120 0v160'/%3E%3C/g%3E%3C/svg%3E"),
    radial-gradient(circle at 70% 30%, #1d8a7a, #0a3d44 60%);
}}
.hero::after {{
  content: "";
  position: absolute;
  inset: auto -10% -20% 35%;
  height: 70%;
  background: radial-gradient(ellipse at center, rgba(217,161,91,.35), transparent 65%);
  pointer-events: none;
  animation: drift 12s ease-in-out infinite alternate;
}}
@keyframes drift {{
  from {{ transform: translate3d(0,0,0) scale(1); }}
  to {{ transform: translate3d(-4%, 3%, 0) scale(1.05); }}
}}
.hero-inner {{
  position: relative;
  z-index: 1;
  max-width: 920px;
}}
.brand {{
  font-family: "Fraunces", serif;
  font-size: clamp(2.6rem, 7vw, 5.2rem);
  font-weight: 700;
  letter-spacing: -.03em;
  line-height: .95;
  margin: 0 0 1rem;
  text-shadow: 0 10px 40px rgba(0,0,0,.25);
}}
.hero h1 {{
  font-family: "Sora", sans-serif;
  font-weight: 500;
  font-size: clamp(1.05rem, 2.2vw, 1.35rem);
  margin: 0 0 .75rem;
  max-width: 38ch;
}}
.hero p.lede {{
  margin: 0 0 1.75rem;
  max-width: 48ch;
  opacity: .92;
  font-size: 1.02rem;
}}
.cta {{
  display: inline-flex;
  gap: .75rem;
  flex-wrap: wrap;
}}
.cta a {{
  text-decoration: none;
  color: var(--sea-deep);
  background: #f7fffc;
  padding: .85rem 1.2rem;
  border-radius: 999px;
  font-weight: 600;
  font-size: .92rem;
  transition: transform .2s ease, box-shadow .2s ease;
  box-shadow: 0 8px 24px rgba(0,0,0,.18);
}}
.cta a.ghost {{
  background: transparent;
  color: #f7fffc;
  border: 1px solid rgba(247,255,252,.45);
  box-shadow: none;
}}
.cta a:hover {{ transform: translateY(-2px); }}
.wrap {{
  width: min(1180px, calc(100% - 2rem));
  margin: 0 auto;
  padding: 2.5rem 0 4rem;
}}
.stats {{
  display: grid;
  grid-template-columns: repeat(4, 1fr);
  gap: .9rem;
  margin: -2.5rem 0 2.5rem;
  position: relative;
  z-index: 2;
}}
.stat {{
  background: var(--card);
  backdrop-filter: blur(8px);
  border: 1px solid var(--line);
  border-radius: 18px;
  padding: 1.1rem 1.15rem;
  box-shadow: 0 10px 30px rgba(18,35,46,.06);
}}
.stat b {{
  display: block;
  font-family: "Fraunces", serif;
  font-size: 2rem;
  line-height: 1;
  color: var(--sea-deep);
}}
.stat span {{ color: var(--muted); font-size: .86rem; }}
.toolbar {{
  display: flex;
  flex-wrap: wrap;
  gap: .6rem;
  align-items: center;
  margin-bottom: 1.5rem;
}}
.toolbar input {{
  flex: 1 1 220px;
  border: 1px solid var(--line);
  border-radius: 999px;
  padding: .75rem 1rem;
  font: inherit;
  background: #fff;
}}
.chip {{
  border: 1px solid var(--line);
  background: #fff;
  border-radius: 999px;
  padding: .55rem .9rem;
  font: inherit;
  cursor: pointer;
  color: var(--muted);
}}
.chip.active, .chip:hover {{
  background: var(--sea-deep);
  color: #fff;
  border-color: var(--sea-deep);
}}
.intro {{
  max-width: 68ch;
  margin-bottom: 2rem;
}}
.intro h2 {{
  font-family: "Fraunces", serif;
  font-size: 1.85rem;
  margin: 0 0 .5rem;
}}
.intro p {{ color: var(--muted); margin: 0; }}
.clade-block {{ margin: 2.25rem 0 2.75rem; }}
.clade-head {{
  display: flex;
  justify-content: space-between;
  gap: 1rem;
  align-items: baseline;
  border-bottom: 2px solid rgba(15,106,106,.25);
  margin-bottom: 1rem;
  padding-bottom: .4rem;
}}
.clade-head h2 {{
  font-family: "Fraunces", serif;
  margin: 0;
  font-size: 1.55rem;
}}
.clade-head p {{ margin: 0; color: var(--muted); font-size: .9rem; }}
.grid {{
  display: grid;
  grid-template-columns: repeat(auto-fill, minmax(280px, 1fr));
  gap: 1rem;
}}
.species {{
  background: var(--card);
  border: 1px solid var(--line);
  border-radius: 18px;
  padding: 1rem 1.05rem 1.15rem;
  box-shadow: 0 8px 24px rgba(18,35,46,.05);
  transition: transform .18s ease, box-shadow .18s ease;
}}
.species:hover {{
  transform: translateY(-3px);
  box-shadow: 0 14px 32px rgba(18,35,46,.1);
}}
.species header {{
  display: flex;
  justify-content: space-between;
  gap: .5rem;
  align-items: center;
  margin-bottom: .55rem;
}}
.pill {{
  font-size: .68rem;
  font-weight: 700;
  letter-spacing: .06em;
  text-transform: uppercase;
  padding: .25rem .55rem;
  border-radius: 999px;
  color: #fff;
}}
.pill.new {{ background: var(--mint); }}
.pill.rescued {{ background: var(--sand); color: #2c2110; }}
.pill.expanded {{ background: var(--sea); }}
.clade-tag {{ font-size: .75rem; color: var(--muted); }}
.species h3 {{
  margin: 0 0 .85rem;
  font-size: 1.12rem;
  font-weight: 600;
}}
.species h3 em {{ font-style: italic; font-family: "Fraunces", serif; font-weight: 600; }}
.metrics {{
  display: grid;
  grid-template-columns: repeat(3, 1fr);
  gap: .45rem;
  margin-bottom: .75rem;
}}
.metric {{
  background: rgba(15,106,106,.06);
  border-radius: 12px;
  padding: .55rem .45rem;
  text-align: center;
}}
.metric b {{ display: block; font-size: 1.05rem; color: var(--sea-deep); }}
.metric span {{ font-size: .68rem; color: var(--muted); }}
.idrow {{
  display: grid;
  grid-template-columns: auto 1fr auto;
  gap: .5rem;
  align-items: center;
  font-size: .78rem;
  color: var(--muted);
  margin-bottom: .65rem;
}}
.idbar {{
  display: block;
  height: 7px;
  border-radius: 999px;
  background: rgba(18,35,46,.08);
  overflow: hidden;
}}
.idbar i {{
  display: block;
  height: 100%;
  background: linear-gradient(90deg, var(--mint), var(--sand));
}}
.meta {{
  display: grid;
  grid-template-columns: 7.2rem 1fr;
  gap: .35rem;
  font-size: .78rem;
  margin-top: .35rem;
}}
.meta span {{ color: var(--muted); }}
.meta code, .meta b {{
  font-family: ui-monospace, SFMono-Regular, Menlo, monospace;
  font-size: .72rem;
  font-weight: 500;
  overflow-wrap: anywhere;
}}
.note {{
  margin: .75rem 0 0;
  padding: .65rem .7rem;
  border-left: 3px solid var(--sand);
  background: rgba(217,161,91,.12);
  border-radius: 0 10px 10px 0;
  font-size: .82rem;
  color: #5a4120;
}}
.expanded-block {{
  margin-top: 3rem;
  padding-top: 1rem;
  border-top: 1px dashed rgba(18,35,46,.2);
}}
.footer {{
  margin-top: 3rem;
  color: var(--muted);
  font-size: .9rem;
}}
.hidden {{ display: none !important; }}
@media (max-width: 800px) {{
  .stats {{ grid-template-columns: repeat(2, 1fr); margin-top: -1.5rem; }}
  .clade-head {{ flex-direction: column; align-items: flex-start; }}
  .hero {{ min-height: 78vh; }}
}}
</style>
</head>
<body>
<header class="hero">
  <div class="hero-inner">
    <p class="brand">Curator Recovery</p>
    <h1>Sixteen genomes where innexin-like loci were newly found or rescued</h1>
    <p class="lede">A genome-direct curator re-probe across discovery species turned sparse automated hits into present, full-length innexin-like models — especially in rotifers, molluscs, spiders and a polychaete.</p>
    <div class="cta">
      <a href="#gallery">Browse species</a>
      <a class="ghost" href="../showcase/index.html">Full thesis showcase</a>
    </div>
  </div>
</header>

<main class="wrap" id="gallery">
  <div class="stats">
    <div class="stat"><b>{n_new}</b><span>new / rescued species</span></div>
    <div class="stat"><b>{n_loci_new}</b><span>present loci in those species</span></div>
    <div class="stat"><b>{n_clades}</b><span>major clades represented</span></div>
    <div class="stat"><b>{n_present}</b><span>total present species ({n_loci} loci)</span></div>
  </div>

  <div class="intro">
    <h2>What “new” means here</h2>
    <p>
      <strong>New</strong> = curator found ≥1 present locus and the species was not previously run, or prior accepted count was 0 (<em>rescued</em>).
      <strong>Expanded</strong> = already had accepted loci, but the curator pass recovered additional / stronger models.
      Present loci are long innexin-like products from the genome-direct probe (not fragmentary-only hits).
    </p>
  </div>

  <div class="toolbar">
    <input id="q" type="search" placeholder="Filter by species name…" aria-label="Filter species">
    <button class="chip active" data-filter="all">All new/rescued</button>
    <button class="chip" data-filter="Rotifera">Rotifera</button>
    <button class="chip" data-filter="Mollusca">Mollusca</button>
    <button class="chip" data-filter="Arthropoda">Arthropoda</button>
    <button class="chip" data-filter="Annelida">Annelida</button>
  </div>

  {''.join(clade_sections)}

  <section class="expanded-block" id="expanded">
    <div class="clade-head">
      <h2>Already known — expanded by curator</h2>
      <p>{len(expanded)} species with prior accepted loci</p>
    </div>
    <div class="grid">
      {expanded_cards}
    </div>
  </section>

  <footer class="footer">
    <p>
      Source: <a href="../gene_curator_probe/innexin_search_summary.md">innexin_search_summary.md</a>
      · Story: <a href="../phylogenetic_story/index.html">phylogenetic story</a>
      · Subfamilies: <a href="../innexin_subfamilies/index.html">innexin subfamilies</a>
    </p>
    <p>Built from curator batch search over 51 genomes · {n_present} with present loci · {n_new} new/rescued ({n_loci_new} loci).</p>
  </footer>
</main>

<script>
const chips = [...document.querySelectorAll('.chip')];
const cards = [...document.querySelectorAll('.species')];
const q = document.getElementById('q');
let filter = 'all';

function apply() {{
  const query = (q.value || '').trim().toLowerCase();
  cards.forEach(card => {{
    const clade = card.dataset.clade;
    const status = card.dataset.status;
    const name = card.dataset.name || '';
    const isNewish = status === 'new' || status === 'rescued';
    const cladeOk = filter === 'all' ? isNewish : (clade === filter && isNewish);
    // expanded cards live in their own section; keep them visible unless searching
    const inExpanded = status === 'expanded';
    const searchOk = !query || name.includes(query);
    let show;
    if (inExpanded) {{
      show = searchOk && (filter === 'all' || clade === filter);
    }} else {{
      show = cladeOk && searchOk;
    }}
    card.classList.toggle('hidden', !show);
  }});
  document.querySelectorAll('.clade-block').forEach(block => {{
    const visible = [...block.querySelectorAll('.species')].some(c => !c.classList.contains('hidden'));
    block.classList.toggle('hidden', !visible);
  }});
}}

chips.forEach(chip => chip.addEventListener('click', () => {{
  chips.forEach(c => c.classList.remove('active'));
  chip.classList.add('active');
  filter = chip.dataset.filter;
  apply();
}}));
q.addEventListener('input', apply);
</script>
</body>
</html>
"""


def main() -> None:
    if not SUMMARY.exists():
        raise SystemExit(f"Missing summary: {SUMMARY}")
    species = _load_species()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    OUT_HTML.write_text(build_html(species), encoding="utf-8")
    new = sum(1 for s in species if s["is_new"])
    print(f"Wrote {OUT_HTML}")
    print(f"Present species={len(species)} new/rescued={new}")


if __name__ == "__main__":
    main()
