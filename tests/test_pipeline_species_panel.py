from __future__ import annotations

import csv
from pathlib import Path

import pipeline.species_panel_builder as panel


class TestSpeciesPanelBuilder:
    def test_build_panel_merges_reference_and_discovery(self, pipeline_sandbox):
        rows = panel.build_panel()
        organisms = {r["organism"] for r in rows}
        assert "Drosophila melanogaster" in organisms
        assert "Test species" in organisms
        roles = {r["panel_role"] for r in rows}
        assert "reference_innexin" in roles
        assert "discovery_innexin" in roles

    def test_build_panel_deduplicates(self, pipeline_sandbox):
        rows = panel.build_panel()
        keys = [(r["family"], r["organism"], r["panel_role"]) for r in rows]
        assert len(keys) == len(set(keys))

    def test_main_writes_csv(self, pipeline_sandbox, tmp_path: Path):
        out = pipeline_sandbox["metadata"] / "panel.csv"
        assert panel.main(["--output", str(out)]) == 0
        rows = list(csv.DictReader(out.open(encoding="utf-8")))
        assert len(rows) >= 2

    def test_load_species_txt_skips_comments(self, tmp_path: Path):
        path = tmp_path / "species.txt"
        path.write_text("# comment\nAlpha species\n\nBeta species\n", encoding="utf-8")
        assert panel.load_species_txt(path) == ["Alpha species", "Beta species"]
