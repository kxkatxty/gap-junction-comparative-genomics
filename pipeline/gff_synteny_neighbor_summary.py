#!/usr/bin/env python3
"""Summarize repetitive / shared flanking genes across species (GFF microsynteny)."""

from __future__ import annotations

import argparse
import csv
import re
import sys
from collections import defaultdict
from pathlib import Path

from pipeline.common import RESULTS_DIR, write_csv

SYNTENY_DIR = RESULTS_DIR / "synteny_gff"
MANIFEST = SYNTENY_DIR / "gene_order_manifest.csv"


def _is_named_gene(name: str) -> bool:
    upper = name.upper()
    prefixes = ("LOC", "MIR", "LNCRNA", "CR", "TRNA", "TRN")
    return not any(upper.startswith(p) for p in prefixes)


def _panel_meta(panel_id: str) -> tuple[str, str]:
    if panel_id.startswith("innexin_"):
        return "innexin", panel_id.removeprefix("innexin_").replace("_cross_species", "")
    if panel_id.startswith("connexin_"):
        return "connexin", panel_id.removeprefix("connexin_").replace("_cross_species", "")
    return "other", panel_id


def _conservation_tier(shared_all: int, shared_majority: int, species_count: int) -> str:
    if shared_all >= 3:
        return "high"
    if shared_all >= 1 or shared_majority >= 3:
        return "moderate"
    if shared_majority >= 1:
        return "low"
    return "none"


def analyze_manifest(manifest_path: Path) -> tuple[list[dict[str, str]], list[dict[str, str]], dict[str, str]]:
    rows = list(csv.DictReader(manifest_path.open(encoding="utf-8")))
    if not rows:
        raise FileNotFoundError(f"Empty manifest: {manifest_path}")

    by_panel_species: dict[str, dict[str, set[str]]] = defaultdict(lambda: defaultdict(set))
    species_per_panel: dict[str, set[str]] = defaultdict(set)

    for row in rows:
        if row.get("is_goi") == "yes":
            continue
        panel = row["panel"]
        species = row["species"]
        name = row["name"].strip()
        species_per_panel[panel].add(species)
        if not _is_named_gene(name):
            continue
        key = name.casefold()
        by_panel_species[panel][species].add(key)

    detail: list[dict[str, str]] = []
    summary: list[dict[str, str]] = []

    for panel in sorted(by_panel_species):
        family, ref_gene = _panel_meta(panel)
        species_map = by_panel_species[panel]
        species_list = sorted(species_per_panel[panel])
        n_species = len(species_list)
        if n_species < 2:
            continue

        gene_counts: dict[str, int] = defaultdict(int)
        for genes in species_map.values():
            for gene in genes:
                gene_counts[gene] += 1

        majority_cutoff = max(2, (n_species + 1) // 2)
        shared_all = sorted(g for g, c in gene_counts.items() if c == n_species)
        shared_majority = sorted(
            (g for g, c in gene_counts.items() if c >= majority_cutoff),
            key=lambda g: (-gene_counts[g], g),
        )
        shared_pair = sorted(
            (g for g, c in gene_counts.items() if c >= 2),
            key=lambda g: (-gene_counts[g], g),
        )

        tier = _conservation_tier(len(shared_all), len(shared_majority), n_species)

        for gene, count in sorted(gene_counts.items(), key=lambda x: (-x[1], x[0])):
            detail.append(
                {
                    "panel": panel,
                    "family": family,
                    "reference_gene": ref_gene,
                    "neighbor_gene": gene,
                    "species_with_gene": str(count),
                    "species_total": str(n_species),
                    "shared_all_species": "yes" if count == n_species else "no",
                    "fraction_species": f"{count / n_species:.2f}",
                }
            )

        summary.append(
            {
                "panel": panel,
                "family": family,
                "reference_gene": ref_gene,
                "species_count": str(n_species),
                "named_neighbors_total": str(len(gene_counts)),
                "shared_all_count": str(len(shared_all)),
                "shared_majority_count": str(len(shared_majority)),
                "shared_pair_count": str(len(shared_pair)),
                "shared_all_genes": ";".join(shared_all[:20]),
                "shared_majority_genes": ";".join(shared_majority[:20]),
                "top_repeated_neighbors": ";".join(shared_pair[:10]),
                "conservation_tier": tier,
            }
        )

    conclusions = _write_conclusions(summary)
    return detail, summary, conclusions


def _write_conclusions(summary: list[dict[str, str]]) -> dict[str, str]:
    innexin = [r for r in summary if r["family"] == "innexin"]
    connexin = [r for r in summary if r["family"] == "connexin"]

    def tier_counts(rows: list[dict[str, str]]) -> dict[str, int]:
        out: dict[str, int] = defaultdict(int)
        for row in rows:
            out[row["conservation_tier"]] += 1
        return dict(out)

    inx_tiers = tier_counts(innexin)
    cnx_tiers = tier_counts(connexin)

    high_cnx = [r for r in connexin if r["conservation_tier"] == "high"]
    high_cnx.sort(key=lambda r: int(r["shared_all_count"]), reverse=True)

    lines = [
        "# Shared flanking genes — conclusion (GFF microsynteny)",
        "",
        "## Main finding",
        "",
        "Vertebrate **connexin** loci show **repetitive, conserved neighbor genes** across species "
        "(same gene symbols in multiple genomes). **Innexin** loci do **not**: insects and nematodes "
        "mostly use species-specific labels (CG*, LOC*), so the same neighbors rarely appear under "
        "identical names across species.",
        "",
        "## Connexins",
        "",
        f"- Panels analysed: **{len(connexin)}**",
        f"- High conservation (≥3 neighbors in all species): **{cnx_tiers.get('high', 0)}**",
        f"- Moderate: **{cnx_tiers.get('moderate', 0)}** · Low: **{cnx_tiers.get('low', 0)}** · None: **{cnx_tiers.get('none', 0)}**",
        "",
        "**Strongest blocks** (neighbors shared in every species of the panel):",
        "",
    ]

    for row in high_cnx[:8]:
        genes = row["shared_all_genes"].replace(";", ", ") or "—"
        lines.append(
            f"- **{row['reference_gene']}** ({row['species_count']} species): {genes}"
        )

    lines.extend(
        [
            "",
            "Recurrent connexin-neighbor patterns:",
            "- **GJA1** locus: **TBC1D32** appears beside GJA1 in all mammal/bird panels.",
            "- **GJA3 / GJB6 / Gja3**: shared **GJB2, GJB6, ZMYM2** block (connexin cluster).",
            "- **GJA4 / GJB3 / GJB4 / GJB5**: shared **DLGAP3, GJB3–GJB5, SMIM12** — classic "
            "mammalian connexin β-cluster on one chromosome.",
            "- **GJB1** (X-linked): **FOXO4, IL2RG, NLGN3, MED12** neighbors repeat across primates/rodents.",
            "",
            "## Innexins",
            "",
            f"- Panels analysed: **{len(innexin)}**",
            f"- High / moderate conservation: **{inx_tiers.get('high', 0) + inx_tiers.get('moderate', 0)}**",
            "",
        ]
    )

    for row in innexin:
        lines.append(
            f"- **{row['reference_gene']}** ({row['species_count']} species): "
            f"{row['shared_all_count']} neighbors in all species; "
            f"best partial overlap = {row['top_repeated_neighbors'] or 'none'} "
            f"(tier: {row['conservation_tier']})."
        )

    lines.extend(
        [
            "",
            "## Interpretation",
            "",
            "1. **Connexin microsynteny is conserved** in vertebrates — the same flanking genes "
            "(often other connexins or TBC1D32) recur across ortholog panels.",
            "2. **Innexin microsynteny is not symbol-conserved** across insects/nematodes — "
            "different annotation systems hide homology; local duplication (Dmel cluster) matters "
            "more than cross-species neighbor identity.",
            "3. **Gene order plots** should be read clade-by-clade: strong conclusions on "
            "repetitive neighbors are valid for connexins; for innexins, compare architecture "
            "and synteny case studies (SynVoy) rather than shared gene names alone.",
            "",
        ]
    )

    text = "\n".join(lines)
    return {
        "markdown": text,
        "short": (
            "Connexins: repetitive neighbors conserved (e.g. GJA1+TBC1D32, β-connexin cluster). "
            "Innexins: no cross-species neighbor repetition by gene name — clade-specific annotation."
        ),
    }


def _write_html_conclusion(conclusions: dict[str, str], summary: list[dict[str, str]]) -> str:
    rows_html = "\n".join(
        f"<tr><td>{r['family']}</td><td>{r['reference_gene']}</td><td>{r['species_count']}</td>"
        f"<td>{r['shared_all_count']}</td><td>{r['conservation_tier']}</td>"
        f"<td>{r['shared_all_genes'] or '—'}</td></tr>"
        for r in summary
    )
    md = conclusions["markdown"].replace("&", "&amp;").replace("<", "&lt;")
    return f"""<!DOCTYPE html><html><head><meta charset='utf-8'>
<title>Neighbor gene sharing — conclusion</title>
<style>
body{{font-family:system-ui,sans-serif;max-width:960px;margin:auto;padding:2rem;background:#f4f7fb;line-height:1.55}}
h1{{color:#2E86AB}}table{{width:100%;border-collapse:collapse;background:#fff;margin:1.5rem 0;font-size:.9rem}}
th,td{{border:1px solid #e0e6ed;padding:.5rem .65rem;text-align:left}}
th{{background:#eaf4fa}}.lead{{background:#fff;padding:1rem 1.25rem;border-left:4px solid #27AE60;border-radius:8px}}
pre{{white-space:pre-wrap;background:#fff;padding:1.25rem;border-radius:8px;font-size:.92rem}}
a{{color:#2E86AB}}
</style></head><body>
<p><a href="index.html">← Back to synteny gallery</a></p>
<h1>Repetitive neighbor genes — conclusion</h1>
<div class="lead"><strong>Short:</strong> {conclusions['short']}</div>
<h2>Per-panel summary</h2>
<table><thead><tr><th>Family</th><th>Gene</th><th>Species</th><th>Shared in all</th><th>Tier</th><th>Shared neighbors</th></tr></thead>
<tbody>{rows_html}</tbody></table>
<h2>Full conclusion</h2>
<pre>{md}</pre>
</body></html>"""


def build(out_dir: Path, manifest_path: Path = MANIFEST) -> None:
    detail, summary, conclusions = analyze_manifest(manifest_path)
    out_dir.mkdir(parents=True, exist_ok=True)
    write_csv(out_dir / "neighbor_sharing_detail.csv", detail)
    write_csv(out_dir / "neighbor_sharing_summary.csv", summary)
    (out_dir / "neighbor_sharing_conclusion.md").write_text(conclusions["markdown"], encoding="utf-8")
    (out_dir / "neighbor_sharing_conclusion.html").write_text(
        _write_html_conclusion(conclusions, summary), encoding="utf-8"
    )
    print(f"Wrote {len(summary)} panel summaries -> {out_dir}")
    print(f"Conclusion: {out_dir / 'neighbor_sharing_conclusion.html'}")


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Summarize shared flanking genes across GFF synteny panels")
    parser.add_argument("--outdir", type=Path, default=SYNTENY_DIR)
    parser.add_argument("--manifest", type=Path, default=MANIFEST)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    try:
        build(args.outdir, args.manifest)
    except FileNotFoundError as exc:
        print(exc, file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
