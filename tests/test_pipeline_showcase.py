"""Tests for scientific showcase builder."""

from __future__ import annotations

from pathlib import Path

from pipeline import scientific_showcase_builder as show
from pipeline.scientific_showcase_builder import SHOWCASE_DIR, build_showcase


def test_build_showcase_creates_plots(tmp_path: Path) -> None:
    out = tmp_path / "showcase"
    records = build_showcase(out)
    existing = [p for _, _, p in records if p.exists()]
    assert len(existing) >= 20
    assert (out / "index.html").exists()
    assert (out / "presenter.html").exists()
    assert (out / "plot_manifest.csv").exists()


def test_showcase_dir_constant() -> None:
    assert "showcase" in str(SHOWCASE_DIR)


class TestShowcaseHelpers:
    def test_infer_clade(self):
        assert show.infer_clade("Homo sapiens") == "mammal"
        assert show.infer_clade("Gallus gallus") == "bird"
        assert show.infer_clade("Drosophila melanogaster") == "insect"
        assert show.infer_clade("Mystery") == "other"

    def test_escape_html(self):
        assert show._escape_html('<a href="x">&') == "&lt;a href=&quot;x&quot;&gt;&amp;"

    def test_clean_plot_title(self):
        title = show._clean_plot_title("exon_map_innexin_compressed")
        assert isinstance(title, str)
        assert title

    def test_read_fasta_length(self, tmp_path: Path):
        fasta = tmp_path / "s.fa"
        fasta.write_text(">a\nACGT\nACGT\n", encoding="utf-8")
        assert show.read_fasta_length(fasta) == 8

    def test_tip_family(self):
        assert show._tip_family("Drosophila_melanogaster|Inx2") == "innexin"
        assert show._tip_family("Homo_sapiens|GJA1") == "connexin"
        assert show._tip_family("some_candidate_locus") == "discovery"
        assert show._tip_family("mystery_label") == "other"
