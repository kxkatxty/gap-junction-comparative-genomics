"""Tests for pannexin pipeline builders (no SynVoy / no network / no IQ-TREE)."""

from __future__ import annotations

import csv
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

from pipeline import limits
from pipeline.common import classify_pannexin_type, write_csv
from pipeline.newick_plot import iter_tips, ladderize, layout, parse_newick
from pipeline import pannexin_annotation_status_builder as annot
from pipeline import pannexin_clade_comparison_builder as clade
from pipeline import pannexin_curator_builder as curator
from pipeline import pannexin_exon_notes_builder as exons
from pipeline import pannexin_insights_builder as insights
from pipeline import pannexin_panel_builder as panel
from pipeline import pannexin_path_builder as path_builder
from pipeline import pannexin_phylogeny_builder as phylo


SIMPLE_TREE = "((A:0.1,B:0.2)90:0.05,C:0.3);"


class TestLimits:
    def test_thread_cap(self, monkeypatch):
        monkeypatch.delenv("PANX_THREADS", raising=False)
        assert limits.n_threads(4) <= limits.MAX_THREADS
        monkeypatch.setenv("PANX_THREADS", "8")
        assert limits.n_threads() == limits.MAX_THREADS
        monkeypatch.setenv("PANX_THREADS", "1")
        assert limits.n_threads() == 1

    def test_heavy_flags_default_off(self, monkeypatch):
        monkeypatch.delenv("PANX_HEAVY_PHYLO", raising=False)
        monkeypatch.delenv("PANX_RUN_MINIPROT", raising=False)
        assert limits.run_heavy_phylogeny() is False
        assert limits.run_genome_curator() is False


class TestNewickPlot:
    def test_parse_and_layout(self):
        root = parse_newick(SIMPLE_TREE)
        ladderize(root)
        layout(root)
        tips = list(iter_tips(root))
        assert len(tips) == 3
        assert {t.name for t in tips} == {"A", "B", "C"}
        assert all(t.x >= 0 for t in tips)

    def test_iqtree_style_bootstrap_labels(self):
        # mirrors IQ-TREE internal bootstrap:label:length pattern
        tree = "((panx|A|PANX1|X:0.1,panx|B|PANX2|Y:0.2)97:0.05,out|Ciona|inx:0.3);"
        root = parse_newick(tree)
        ladderize(root)
        layout(root)
        tips = list(iter_tips(root))
        assert len(tips) == 3
        assert any(t.name and t.name.startswith("out|") for t in tips)


class TestPanelHelpers:
    def test_false_positive_filter(self):
        assert panel.FALSE_POSITIVE_GENES.search("pknox2")
        assert panel.FALSE_POSITIVE_GENES.search("LRRC8A")
        assert not panel.FALSE_POSITIVE_GENES.search("PANX1")

    def test_needs_type_rescue(self):
        assert panel._needs_type_rescue({"gene": "PANX", "panx_type": "other/unknown"})
        assert panel._needs_type_rescue({"gene": "PANX", "panx_type": "PANX2"})
        assert not panel._needs_type_rescue({"gene": "LOC1", "panx_type": "PANX3"})
        assert not panel._needs_type_rescue({"gene": "Panx1", "panx_type": "PANX1"})

    def test_dedupe_keeps_longest(self):
        rows = [
            {
                "organism_key": "Homo_sapiens",
                "organism": "Homo sapiens",
                "gene": "PANX1",
                "accession": "A",
                "panx_type": "PANX1",
                "clade": "Mammal",
                "length": "100",
                "seq_id": "a",
                "seq": "A" * 100,
            },
            {
                "organism_key": "Homo_sapiens",
                "organism": "Homo sapiens",
                "gene": "PANX1",
                "accession": "B",
                "panx_type": "PANX1",
                "clade": "Mammal",
                "length": "400",
                "seq_id": "b",
                "seq": "A" * 400,
            },
        ]
        out = panel.dedupe_longest(rows)
        assert len(out) == 1
        assert out[0]["accession"] == "B"

    def test_write_fasta_and_summarize(self, tmp_path, monkeypatch):
        monkeypatch.setattr(panel, "OUT", tmp_path)
        monkeypatch.setattr(panel, "CLADE_SUMMARY", tmp_path / "clade.csv")
        monkeypatch.setattr(panel, "TYPE_SUMMARY", tmp_path / "type.csv")
        rows = [
            {
                "seq_id": "panx|Mus|Panx1|Q",
                "organism": "Mus musculus",
                "organism_key": "Mus_musculus",
                "gene": "Panx1",
                "accession": "Q",
                "panx_type": "PANX1",
                "clade": "Mammal",
                "length": "426",
                "cys_count": "10",
                "source_path": "x",
                "seq": "M" * 426,
            }
        ]
        fa = tmp_path / "all.fasta"
        panel.write_fasta(rows, fa)
        assert fa.read_text().startswith(">panx|Mus|Panx1|Q")
        panel.summarize(rows)
        assert (tmp_path / "clade.csv").exists()
        assert (tmp_path / "type.csv").exists()


class TestPhylogenyHelpers:
    def test_tip_type_and_label(self):
        meta = {
            "panx|Homo_sapiens|PANX1|Q96RD7": {
                "organism": "Homo sapiens",
                "gene": "PANX1",
                "panx_type": "PANX1",
            }
        }
        assert phylo._tip_type("panx|Homo_sapiens|PANX1|Q96RD7", meta) == "PANX1"
        assert phylo._tip_type("out|Ciona|inx1", meta) == "outgroup_inx"
        lab = phylo._short_label("panx|Homo_sapiens|PANX1|Q96RD7", meta)
        assert lab.startswith("1:")
        assert "PANX1" in lab or "Homo" in lab

    def test_plot_iqtree_tree(self, tmp_path, monkeypatch):
        tree = tmp_path / "t.treefile"
        tree.write_text("((panx|A|PANX1|X:0.1,panx|B|PANX2|Y:0.2)50:0.05,out|Ciona|inx:0.3);", encoding="utf-8")
        monkeypatch.setattr(phylo, "TREE", tree)
        meta = [
            {"seq_id": "panx|A|PANX1|X", "organism": "A a", "gene": "PANX1", "panx_type": "PANX1"},
            {"seq_id": "panx|B|PANX2|Y", "organism": "B b", "gene": "PANX2", "panx_type": "PANX2"},
        ]
        out = tmp_path / "tree.png"
        phylo.plot_iqtree_tree(meta, out)
        assert out.exists() and out.stat().st_size > 1000

    def test_write_html_no_env_jargon(self, tmp_path, monkeypatch):
        monkeypatch.setattr(phylo, "OUT", tmp_path)
        monkeypatch.setattr(phylo, "TREE", tmp_path / "missing.treefile")
        phylo.write_html(10, 2, "LG+G4")
        html = (tmp_path / "index.html").read_text(encoding="utf-8")
        assert "IQ-TREE" in html
        assert "PANX_HEAVY" not in html
        assert "rendered" in html.lower() or "phylogram" in html.lower() or "treefile" in html.lower()


class TestAnnotationStatus:
    def test_false_positive_hit_pattern(self):
        assert annot.FALSE_POSITIVE_HIT.search("pknox2")
        assert annot.FALSE_POSITIVE_HIT.search("GRB10-interacting GYF protein 2")
        assert not annot.FALSE_POSITIVE_HIT.search("Pannexin-1")

    def test_enrich_filters_pknox(self, tmp_path, monkeypatch):
        meta = tmp_path / "meta"
        meta.mkdir()
        status = meta / "pannexin_annotation_status.csv"
        write_csv(
            status,
            [
                {
                    "organism": "Danio rerio",
                    "family": "pannexin",
                    "final_category": "annotated_family_member_found",
                    "confidence": "high",
                    "top_matching_term": "pannexin",
                    "top_gene_symbol": "pknox2",
                    "top_protein_name": "Pannexin",
                    "database_source": "uniprot",
                }
            ],
        )
        cfg = meta / "species_config.csv"
        write_csv(
            cfg,
            [{"family": "pannexin", "organism": "Danio rerio", "reviewed_only": "false", "clade_hint": "Ray-finned fish"}],
        )
        panel_dir = tmp_path / "panel"
        panel_dir.mkdir()
        write_csv(
            panel_dir / "sequence_metadata.csv",
            [
                {
                    "seq_id": "x",
                    "organism": "Danio rerio",
                    "panx_type": "PANX1",
                    "clade": "Ray-finned fish",
                },
                {
                    "seq_id": "y",
                    "organism": "Danio rerio",
                    "panx_type": "PANX2",
                    "clade": "Ray-finned fish",
                },
            ],
        )
        out = tmp_path / "annot"
        out.mkdir()
        monkeypatch.setattr(annot, "METADATA_DIR", meta)
        monkeypatch.setattr(annot, "SPECIES_CONFIG", cfg)
        monkeypatch.setattr(annot, "STATUS_CSV", status)
        monkeypatch.setattr(annot, "OUT", out)
        monkeypatch.setattr(annot, "CLADE_SUMMARY", out / "status_by_clade.csv")
        monkeypatch.setattr(annot, "PANEL_META", panel_dir / "sequence_metadata.csv")
        rows = annot.enrich_and_summarize()
        assert len(rows) == 1
        assert rows[0]["organism"] == "Danio rerio"
        assert "filtered" in rows[0]["qc_note"].lower()
        assert rows[0]["top_matching_term"] == "panel_curated_pannexin"
        assert rows[0]["n_panel_proteins"] == "2"


class TestCladeComparison:
    def test_load_rows_from_meta(self, tmp_path, monkeypatch):
        meta = tmp_path / "sequence_metadata.csv"
        write_csv(
            meta,
            [
                {
                    "seq_id": "a",
                    "organism": "Homo sapiens",
                    "gene": "PANX1",
                    "accession": "Q",
                    "panx_type": "PANX1",
                    "clade": "Mammal",
                    "length": "426",
                    "cys_count": "11",
                }
            ],
        )
        monkeypatch.setattr(clade, "PANEL_META", meta)
        rows = clade.load_rows()
        assert len(rows) == 1
        assert rows[0]["length"] == 426
        assert rows[0]["bottleneck_clade"] == "no"

    def test_write_tables(self, tmp_path, monkeypatch):
        monkeypatch.setattr(clade, "OUT", tmp_path)
        rows = [
            {
                "seq_id": "a",
                "organism": "Homo sapiens",
                "organism_key": "Homo_sapiens",
                "gene": "PANX1",
                "accession": "Q",
                "panx_type": "PANX1",
                "clade": "Mammal",
                "length": 426,
                "cys_count": 11,
                "source": "reference_panel",
                "bottleneck_clade": "no",
            },
            {
                "seq_id": "b",
                "organism": "Petromyzon marinus",
                "organism_key": "Petromyzon_marinus",
                "gene": "LOC1",
                "accession": "Z",
                "panx_type": "PANX3",
                "clade": "Jawless vertebrate",
                "length": 450,
                "cys_count": 10,
                "source": "reference_panel",
                "bottleneck_clade": "yes",
            },
        ]
        clade.write_tables(rows)
        assert (tmp_path / "species_copy_number.csv").exists()
        assert (tmp_path / "clade_by_type_matrix.csv").exists()
        findings = clade.key_findings(rows)
        assert findings
        assert not any("no SynVoy" in f for f in findings)


class TestCuratorPathInsightsExons:
    def test_curator_builds(self, tmp_path, monkeypatch):
        monkeypatch.setattr(curator, "OUT", tmp_path)
        monkeypatch.setattr(curator, "STATUS_CSV", tmp_path / "annotation_gap_status.csv")
        curator.build()
        html = (tmp_path / "index.html").read_text(encoding="utf-8")
        assert "Branchiostoma" in html
        assert "other/unknown" not in html  # outdated wording removed
        assert "Callorhinchus" in html

    def test_path_html_has_chapters(self, tmp_path, monkeypatch):
        monkeypatch.setattr(path_builder, "OUT", tmp_path)
        monkeypatch.setattr(path_builder, "FIGS", tmp_path / "figures")
        (tmp_path / "figures").mkdir()
        # skip figure copy warnings by touching expected names
        for name in ("05_bottleneck_copies.png", "06_clade_type_heatmap.png", "07_type_mix.png"):
            (tmp_path / "figures" / name).write_bytes(b"")
        path_builder.write_html()
        html = (tmp_path / "index.html").read_text(encoding="utf-8")
        assert "Across clades" in html
        assert "(done)" not in html
        assert "species_config.csv" not in html

    def test_insights_build_with_missing_figs(self, tmp_path, monkeypatch):
        monkeypatch.setattr(insights, "OUT", tmp_path)
        monkeypatch.setattr(insights, "FIGS", tmp_path / "figures")
        monkeypatch.setattr(insights, "RESULTS_DIR", tmp_path)
        # empty panel meta
        panel_dir = tmp_path / "pannexin_panel"
        panel_dir.mkdir()
        write_csv(panel_dir / "sequence_metadata.csv", [])
        (tmp_path / "panx_vs_inx").mkdir()
        write_csv(tmp_path / "panx_vs_inx" / "summary.csv", [{"pct_with_hit": "90%"}])
        (tmp_path / "pannexin_clade_comparison").mkdir()
        write_csv(tmp_path / "pannexin_clade_comparison" / "species_copy_number.csv", [])
        insights.build()
        html = (tmp_path / "index.html").read_text(encoding="utf-8")
        assert "Pannexin insights" in html
        assert "Literature vs this panel" in html

    def test_exon_html(self, tmp_path, monkeypatch):
        monkeypatch.setattr(exons, "OUT", tmp_path)
        monkeypatch.setattr(exons, "FIGS", tmp_path / "figures")
        (tmp_path / "figures").mkdir()
        rows = [
            {
                "organism": "Homo sapiens",
                "ensembl_species": "human",
                "symbol": "PANX1",
                "panx_type": "PANX1",
                "gene_id": "ENSG1",
                "n_coding_exons": "5",
                "status": "ok",
            }
        ]
        exons.plot_exons(rows, tmp_path / "figures" / "01_exon_counts.png")
        exons.write_html(rows)
        html = (tmp_path / "index.html").read_text(encoding="utf-8")
        assert "PANX1" in html
        assert "Ensembl" in html


class TestClassifySmoke:
    def test_pannexin_1a_style(self):
        assert classify_pannexin_type("Pannexin-1") == "PANX1"
        assert classify_pannexin_type("pannexin 2-like") == "PANX2"
