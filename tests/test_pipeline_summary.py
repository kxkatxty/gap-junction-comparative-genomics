from __future__ import annotations

import csv

import pipeline.summary_table_builder as summary


class TestSummaryTableBuilder:
    def test_build_master_summary_joins_sources(self, pipeline_sandbox):
        rows = summary.build_master_summary()
        test_row = next(r for r in rows if r["organism"] == "Test species")
        assert test_row["family"] == "innexin"
        assert test_row["annotation_category"] == "no_related_entry_found"
        assert test_row["discovery_accepted"] == "1"
        assert test_row["discovery_high_confidence"] == "1"

        ref_row = next(r for r in rows if r["organism"] == "Drosophila melanogaster")
        assert ref_row["assembly_level"] == "Chromosome"
        assert ref_row["has_gff_annotation"] == "yes"
        assert ref_row["annotation_category"] == "annotated_family_member_found"

    def test_main_writes_csv(self, pipeline_sandbox):
        out = pipeline_sandbox["metadata"] / "master.csv"
        assert summary.main(["--output", str(out)]) == 0
        rows = list(csv.DictReader(out.open(encoding="utf-8")))
        assert any(r["organism"] == "Test species" for r in rows)
