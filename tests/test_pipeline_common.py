from __future__ import annotations

import csv
import json
from pathlib import Path

import pipeline.common as common


class TestCommonHelpers:
    def test_read_csv_rows_missing_file(self, tmp_path: Path):
        assert common.read_csv_rows(tmp_path / "missing.csv") == []

    def test_read_csv_rows_empty_file(self, tmp_path: Path):
        path = tmp_path / "empty.csv"
        path.write_text("", encoding="utf-8")
        assert common.read_csv_rows(path) == []

    def test_write_csv_roundtrip(self, tmp_path: Path):
        out = tmp_path / "out.csv"
        rows = [{"a": "1", "b": "2"}, {"a": "3", "b": "4"}]
        common.write_csv(out, rows)
        assert common.read_csv_rows(out) == rows

    def test_write_csv_empty_writes_empty_file(self, tmp_path: Path):
        out = tmp_path / "empty.csv"
        common.write_csv(out, [])
        assert out.read_text(encoding="utf-8") == ""

    def test_infer_family_from_path(self):
        assert common.infer_family_from_path(Path("project/results/innexin_discovery/x")) == "innexin"
        assert common.infer_family_from_path(Path("project/results/connexin_evolutionary_discovery/x")) == "connexin"
        assert common.infer_family_from_path(Path("project/results/other")) == "unknown"

    def test_iter_discovery_summaries(self, pipeline_sandbox):
        summaries = common.iter_discovery_summaries()
        assert len(summaries) == 1
        path, data, family = summaries[0]
        assert family == "innexin"
        assert data["organism"] == "Test species"
        assert path.name == "discovery_summary.json"
