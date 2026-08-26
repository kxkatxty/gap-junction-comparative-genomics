#!/usr/bin/env python3
"""GFF-based microsynteny for gap-junction genes (no SynVoy).

Builds gene-order plots around each annotated locus in gene_summary.csv using
NCBI GFF annotations. Cross-species panels group orthologs by reference symbol.
"""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from pipeline.common import ANNOTATIONS_DIR, RESULTS_DIR, write_csv
from pipeline.foxp2_style_synteny import (
    LocusSpec,
    build_panel,
    default_panels,
    load_row,
    plot_panel,
)

OUT_DIR = RESULTS_DIR / "synteny_gff"
GENE_SUMMARY = RESULTS_DIR / "exon_structures" / "gene_summary.csv"


def _annotation_gff(family: str, organism: str) -> Path | None:
    folder = organism.replace(" ", "_")
    root = ANNOTATIONS_DIR / family / folder
    if not root.is_dir():
        return None
    matches = sorted(root.glob("*.gff"))
    return matches[0] if matches else None


def _sanitize_id(text: str) -> str:
    slug = re.sub(r"[^A-Za-z0-9]+", "_", text.strip()).strip("_")
    return slug[:80] or "panel"


def locus_from_row(row: dict[str, str]) -> LocusSpec | None:
    gff = _annotation_gff(row["family"], row["organism"])
    if gff is None:
        return None
    goi_ids = tuple(
        x
        for x in (row["gene_id"], row["gene_symbol"], row["reference_gene_symbol"])
        if x and x != "-"
    )
    if not goi_ids:
        return None
    label = f"{row['organism']} — {row['reference_gene_symbol']}"
    return LocusSpec(
        species_label=label,
        gff_path=gff,
        chrom=row["seqid"],
        goi_gene_ids=goi_ids,
        clade=row["family"],
    )


def load_gene_summary(path: Path) -> list[dict[str, str]]:
    if not path.exists():
        return []
    with path.open(encoding="utf-8", newline="") as handle:
        return list(csv.DictReader(handle))


def panels_from_gene_summary(
    rows: list[dict[str, str]],
    *,
    min_species: int = 2,
) -> list[tuple[str, str, list[LocusSpec]]]:
    """Cross-species panels keyed by family + reference_gene_symbol."""
    grouped: dict[tuple[str, str], list[LocusSpec]] = defaultdict(list)
    for row in rows:
        spec = locus_from_row(row)
        if spec is None:
            continue
        _, genes = load_row(spec)
        if not genes or not any(g.is_goi for g in genes):
            continue
        key = (row["family"], row["reference_gene_symbol"])
        grouped[key].append(spec)

    panels: list[tuple[str, str, list[LocusSpec]]] = []
    for (family, ref_symbol), specs in sorted(grouped.items()):
        if len(specs) < min_species:
            continue
        panel_id = _sanitize_id(f"{family}_{ref_symbol}_cross_species")
        title = f"{ref_symbol} ({family}) — microsynteny across {len(specs)} species"
        panels.append((panel_id, title, specs))
    return panels


def single_locus_panels(rows: list[dict[str, str]]) -> list[tuple[str, str, list[LocusSpec]]]:
    panels: list[tuple[str, str, list[LocusSpec]]] = []
    for row in rows:
        spec = locus_from_row(row)
        if spec is None:
            continue
        _, genes = load_row(spec)
        if not genes or not any(g.is_goi for g in genes):
            continue
        panel_id = _sanitize_id(
            f"{row['family']}_{row['organism']}_{row['reference_gene_symbol']}_locus"
        )
        title = f"{row['reference_gene_symbol']} — {row['organism']} ({row['family']})"
        panels.append((panel_id, title, [spec]))
    return panels


def build_all(
    out_dir: Path,
    *,
    mode: str = "all",
    gene_summary_path: Path = GENE_SUMMARY,
) -> list[dict[str, str]]:
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = load_gene_summary(gene_summary_path)
    if not rows:
        raise FileNotFoundError(f"No gene summary at {gene_summary_path}")

    thesis_ids = {p[0] for p in default_panels()}
    panels: list[tuple[str, str, list[LocusSpec]]] = []
    if mode in ("all", "cross"):
        panels.extend(panels_from_gene_summary(rows))
    if mode in ("all", "single"):
        panels.extend(single_locus_panels(rows))
    if mode in ("all", "thesis"):
        panels.extend(default_panels())

    # Deduplicate panel ids
    seen: set[str] = set()
    unique_panels: list[tuple[str, str, list[LocusSpec]]] = []
    for panel_id, title, specs in panels:
        if panel_id in seen:
            continue
        seen.add(panel_id)
        unique_panels.append((panel_id, title, specs))

    results: list[dict[str, str]] = []
    manifest: list[dict[str, str]] = []
    cross_dir = out_dir / "cross_species"
    single_dir = out_dir / "single_locus"
    thesis_dir = out_dir / "case_studies"

    for panel_id, title, specs in unique_panels:
        if panel_id in thesis_ids:
            target = thesis_dir
        elif panel_id.endswith("_cross_species"):
            target = cross_dir
        else:
            target = single_dir
        try:
            meta = build_panel(panel_id, title, specs, target)
            meta["mode"] = target.name
            meta["family"] = specs[0].clade if specs else _panel_family(panel_id)
            results.append(meta)
            print(f"OK  {panel_id} ({target.name})")
            for spec, genes in (load_row(s) for s in specs):
                for g in genes:
                    manifest.append(
                        {
                            "panel": panel_id,
                            "mode": target.name,
                            "species": spec.species_label,
                            "family": spec.clade,
                            "gene_id": g.gene_id,
                            "name": g.name,
                            "chrom": spec.chrom,
                            "start": str(g.start),
                            "end": str(g.end),
                            "strand": g.strand,
                            "is_goi": "yes" if g.is_goi else "no",
                        }
                    )
        except (ValueError, OSError) as exc:
            print(f"SKIP {panel_id}: {exc}", file=sys.stderr)

    write_csv(out_dir / "gene_order_manifest.csv", manifest)
    _write_index(out_dir, results)
    _write_summary(out_dir, rows, results)
    return results


def _write_summary(out_dir: Path, rows: list[dict[str, str]], results: list[dict[str, str]]) -> None:
    summary = [
        {"metric": "genes_in_summary", "value": str(len(rows))},
        {"metric": "panels_built", "value": str(len(results))},
        {
            "metric": "cross_species_panels",
            "value": str(sum(1 for r in results if r.get("mode") == "cross_species")),
        },
        {
            "metric": "single_locus_panels",
            "value": str(sum(1 for r in results if r.get("mode") == "single_locus")),
        },
        {
            "metric": "case_study_panels",
            "value": str(sum(1 for r in results if r.get("mode") == "case_studies")),
        },
    ]
    write_csv(out_dir / "build_summary.csv", summary)


def _panel_family(panel_id: str) -> str:
    if panel_id.startswith("innexin_"):
        return "innexin"
    if panel_id.startswith("connexin_"):
        return "connexin"
    return "other"


def _collect_results_from_disk(out_dir: Path) -> list[dict[str, str]]:
    """Rebuild result rows from existing PNG/HTML on disk (fast index refresh)."""
    results: list[dict[str, str]] = []
    for mode in ("cross_species", "single_locus", "case_studies"):
        folder = out_dir / mode
        if not folder.is_dir():
            continue
        for png in sorted(folder.glob("*.png")):
            html = png.with_suffix(".html")
            results.append(
                {
                    "panel": png.stem,
                    "png": str(png),
                    "html": str(html) if html.exists() else str(png),
                    "mode": mode,
                    "family": _panel_family(png.stem),
                }
            )
    return results


def _write_index(out_dir: Path, results: list[dict[str, str]]) -> None:
    by_mode: dict[str, list[tuple[str, str, str, str]]] = defaultdict(list)
    for row in results:
        mode = row.get("mode", "single_locus")
        rel_html = f"{mode}/{Path(row['html']).name}"
        rel_png = f"{mode}/{Path(row['png']).name}"
        family = row.get("family") or _panel_family(row["panel"])
        by_mode[mode].append((row["panel"], rel_html, rel_png, family))

    def cards(items: list[tuple[str, str, str, str]]) -> str:
        return "\n".join(
            f"<div class='card'><a href='{html}'><img src='{png}' alt='{pid}'><p>{pid.replace('_', ' ')}</p></a></div>"
            for pid, html, png, _family in sorted(items)
        )

    cross = by_mode.get("cross_species", [])
    innexin_cross = [item for item in cross if item[3] == "innexin"]
    connexin_cross = [item for item in cross if item[3] == "connexin"]
    other_cross = [item for item in cross if item[3] not in ("innexin", "connexin")]

    index = out_dir / "index.html"
    index.write_text(
        f"""<!DOCTYPE html><html><head><meta charset='utf-8'>
<title>GFF microsynteny — gap junction genes</title>
<style>
body{{font-family:system-ui,sans-serif;max-width:1100px;margin:auto;padding:2rem;background:#f4f7fb;color:#1a1d26}}
h1{{color:#2E86AB}}h2{{margin-top:2.5rem;color:#34495e}}
.nav{{display:flex;gap:1rem;flex-wrap:wrap;margin:1rem 0 2rem}}
.nav a{{color:#2E86AB;font-weight:600;text-decoration:none;padding:.35rem .75rem;background:#fff;border-radius:6px;border:1px solid #d5e8f0}}
.nav a:hover{{background:#eaf4fa}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:1.5rem}}
.card img{{width:100%;border-radius:8px;box-shadow:0 2px 12px rgba(0,0,0,.08)}}
.card p{{text-align:center;font-weight:600;font-size:.9rem}}
a{{text-decoration:none;color:inherit}}
.note{{background:#fff;padding:1rem 1.25rem;border-radius:8px;border-left:4px solid #2E86AB}}
.count{{color:#5d6d7e;font-size:.95rem;margin:-1rem 0 1rem}}
</style></head><body>
<h1>GFF microsynteny (no SynVoy)</h1>
<div class="note">
<p>Gene-order context from NCBI GFF annotations around each gap-junction locus in
<code>gene_summary.csv</code>. Red = gene of interest; blue = flanking annotated genes.</p>
<p><strong>{len(results)}</strong> panels · innexin: <strong>{len(innexin_cross)}</strong> · connexin: <strong>{len(connexin_cross)}</strong></p>
</div>
<nav class="nav">
<a href="#innexin">Innexins</a>
<a href="#connexin">Connexins</a>
<a href="#conclusion">Conclusion</a>
<a href="#case_studies">Case studies</a>
<a href="#single_locus">Single locus</a>
<a href="neighbor_sharing_conclusion.html">Neighbor sharing report</a>
</nav>
<div id="conclusion" class="note" style="margin-bottom:2rem;border-left-color:#27AE60">
<p><strong>Repetitive neighbors:</strong> Vertebrate connexins share the same flanking gene symbols across species
(e.g. GJA1 + TBC1D32; β-connexin cluster GJB3–GJB5). Innexins show almost no cross-species neighbor repetition by name.
<a href="neighbor_sharing_conclusion.html">Full conclusion →</a></p>
</div>
<h2 id="innexin">Innexins — cross-species</h2>
<p class="count">{len(innexin_cross)} panel(s)</p>
<div class="grid">{cards(innexin_cross) or "<p>No innexin panels yet.</p>"}</div>
<h2 id="connexin">Connexins — cross-species</h2>
<p class="count">{len(connexin_cross)} panel(s)</p>
<div class="grid">{cards(connexin_cross) or "<p>No connexin panels yet.</p>"}</div>
{f'<h2>Other cross-species panels</h2><div class="grid">{cards(other_cross)}</div>' if other_cross else ''}
<h2 id="case_studies">Case studies (Dmel cluster, Inx2 cross-clade)</h2>
<div class="grid">{cards(by_mode.get('case_studies', []))}</div>
<h2 id="single_locus">Single-locus panels (all annotated genes)</h2>
<div class="grid">{cards(by_mode.get('single_locus', []))}</div>
</body></html>""",
        encoding="utf-8",
    )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="GFF microsynteny for gap-junction genes (no SynVoy)")
    parser.add_argument("--outdir", type=Path, default=OUT_DIR)
    parser.add_argument(
        "--mode",
        choices=("all", "cross", "single", "thesis"),
        default="all",
        help="cross = ortholog panels; single = per-gene loci; thesis = curated case studies",
    )
    parser.add_argument("--gene-summary", type=Path, default=GENE_SUMMARY)
    parser.add_argument(
        "--index-only",
        action="store_true",
        help="Rebuild index.html from existing PNGs (no plot regeneration)",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    if args.index_only:
        results = _collect_results_from_disk(args.outdir)
        if not results:
            print(f"No panels found under {args.outdir}", file=sys.stderr)
            return 1
        _write_index(args.outdir, results)
        print(f"Index updated: {args.outdir / 'index.html'} ({len(results)} panels)")
        return 0
    try:
        build_all(args.outdir, mode=args.mode, gene_summary_path=args.gene_summary)
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
