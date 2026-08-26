"""Tests for FOXP2-style synteny plots."""

from __future__ import annotations

from pathlib import Path

from pipeline.foxp2_style_synteny import build_panel, default_panels, parse_genes_on_chrom
from pipeline.common import ANNOTATIONS_DIR


def test_parse_dmel_genes_on_chrom() -> None:
    gff = ANNOTATIONS_DIR / "innexin" / "Drosophila_melanogaster"
    matches = list(gff.glob("*.gff"))
    assert matches
    genes = parse_genes_on_chrom(matches[0], "NC_004354.4")
    assert any(g.name == "Inx2" for g in genes)


def test_build_proximal_panel(tmp_path: Path) -> None:
    panels = {p[0]: p for p in default_panels()}
    panel_id, title, specs = panels["dmel_innexin_cluster_proximal"]
    meta = build_panel(panel_id, title, specs, tmp_path)
    assert Path(meta["png"]).exists()
    assert Path(meta["html"]).exists()
