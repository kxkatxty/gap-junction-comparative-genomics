from __future__ import annotations

from pathlib import Path

import collect_unannotated_connexin_by_clade as cuc
from classify_family_annotation import CATEGORY_ANNOTATED


class TestConnexinNotAnnotated:
    def test_annotated_is_false(self):
        assert cuc.connexin_not_annotated(CATEGORY_ANNOTATED) is False

    def test_missing_is_true(self):
        assert cuc.connexin_not_annotated("no_related_entry_found") is True


class TestLoadExcludeOrganisms:
    def test_load_from_txt(self, tmp_path: Path):
        txt = tmp_path / "exclude.txt"
        txt.write_text("# skip\nAlpha species\n\nBeta species\n", encoding="utf-8")
        excluded = cuc.load_exclude_organisms([txt])
        assert "alphaspecies" in excluded
        assert "betaspecies" in excluded

    def test_load_from_csv(self, tmp_path: Path):
        csv_path = tmp_path / "exclude.csv"
        csv_path.write_text("organism\nGamma species\n", encoding="utf-8")
        excluded = cuc.load_exclude_organisms([csv_path])
        assert "gammaspecies" in excluded

    def test_missing_file_ignored(self, tmp_path: Path):
        assert cuc.load_exclude_organisms([tmp_path / "missing.txt"]) == set()


class TestRankOrganisms:
    def test_skips_model_species_and_low_assembly(self):
        assemblies = {
            "Homo sapiens": {"AssemblyStatus": "Chromosome"},
            "Low quality species": {"AssemblyStatus": "Contig"},
            "Good species": {"AssemblyStatus": "Chromosome"},
        }
        ranked = cuc.rank_organisms(
            assemblies,
            model_species={cuc.normalize_key("Homo sapiens")},
            min_assembly_level="Scaffold",
        )
        organisms = [item[0] for item in ranked]
        assert "Homo sapiens" not in organisms
        assert "Low quality species" not in organisms
        assert organisms == ["Good species"]
