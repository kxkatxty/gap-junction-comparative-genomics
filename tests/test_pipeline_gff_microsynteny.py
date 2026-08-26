"""Tests for GFF microsynteny (no SynVoy)."""

from __future__ import annotations

from pathlib import Path

from pipeline.gff_microsynteny import (
    load_gene_summary,
    locus_from_row,
    panels_from_gene_summary,
    single_locus_panels,
)
from pipeline.foxp2_style_synteny import build_panel, load_row


def test_locus_from_dmel_inx2() -> None:
    rows = load_gene_summary(Path("project/results/exon_structures/gene_summary.csv"))
    row = next(r for r in rows if r["reference_gene_symbol"] == "Inx2" and "melanogaster" in r["organism"])
    spec = locus_from_row(row)
    assert spec is not None
    _, genes = load_row(spec)
    assert any(g.is_goi for g in genes)


def test_cross_species_gja1_panel(tmp_path: Path) -> None:
    rows = load_gene_summary(Path("project/results/exon_structures/gene_summary.csv"))
    panels = panels_from_gene_summary(rows)
    gja1 = next(p for p in panels if "GJA1" in p[0])
    panel_id, title, specs = gja1
    assert len(specs) >= 2
    meta = build_panel(panel_id, title, specs, tmp_path)
    assert Path(meta["png"]).exists()


def test_single_locus_count() -> None:
    rows = load_gene_summary(Path("project/results/exon_structures/gene_summary.csv"))
    singles = single_locus_panels(rows)
    assert len(singles) >= 50
