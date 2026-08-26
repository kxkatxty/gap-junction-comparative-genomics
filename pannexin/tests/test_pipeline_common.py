from __future__ import annotations

from pipeline.common import classify_pannexin_type, write_csv, read_csv_rows


class TestPannexinCommon:
    def test_classify_types(self):
        assert classify_pannexin_type("PANX1_HUMAN") == "PANX1"
        assert classify_pannexin_type("Panx2") == "PANX2"
        assert classify_pannexin_type("PX3") == "PANX3"
        assert classify_pannexin_type("unknown") == "other/unknown"

    def test_csv_roundtrip(self, tmp_path):
        path = tmp_path / "t.csv"
        write_csv(path, [{"a": "1", "b": "2"}])
        assert read_csv_rows(path) == [{"a": "1", "b": "2"}]


class TestPannexinImports:
    def test_all_builders_import(self):
        import pipeline.pannexin_panel_builder  # noqa: F401
        import pipeline.pannexin_clade_comparison_builder  # noqa: F401
        import pipeline.pannexin_similarity_builder  # noqa: F401
        import pipeline.panx_vs_inx_builder  # noqa: F401
        import pipeline.pannexin_phylogeny_builder  # noqa: F401
        import pipeline.pannexin_annotation_status_builder  # noqa: F401
        import pipeline.pannexin_curator_builder  # noqa: F401
        import pipeline.pannexin_path_builder  # noqa: F401
        import pipeline.pannexin_insights_builder  # noqa: F401
        import pipeline.pannexin_exon_notes_builder  # noqa: F401
        import pipeline.newick_plot  # noqa: F401
        import pipeline.limits  # noqa: F401
