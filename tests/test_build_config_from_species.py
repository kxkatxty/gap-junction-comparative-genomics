from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import patch

import build_config_from_species as bcfg


class TestLoadSpeciesNames:
    def test_inline_and_file_deduplicated(self, tmp_path: Path):
        species_file = tmp_path / "species.txt"
        species_file.write_text("Danio rerio\n# comment\nMus musculus\n", encoding="utf-8")
        names = bcfg.load_species_names([species_file], ["Homo sapiens", "Danio rerio"])
        assert names == ["Homo sapiens", "Danio rerio", "Mus musculus"]

    def test_inline_semicolon_separated(self):
        names = bcfg.load_species_names([], ["Homo sapiens; Mus musculus"])
        assert len(names) == 2


class TestTaxonomy:
    def test_lineage_has_vertebrata_by_taxon_id(self):
        lineage = [{"taxonId": 7742, "scientificName": "Vertebrata"}]
        assert bcfg._lineage_has_vertebrata(lineage)

    def test_lineage_has_vertebrata_by_name(self):
        lineage = [{"scientificName": "Vertebrata"}]
        assert bcfg._lineage_has_vertebrata(lineage)

    def test_lineage_no_vertebrata(self):
        assert bcfg._lineage_has_vertebrata([{"scientificName": "Metazoa"}]) is False

    def test_pick_taxonomy_match_exact(self):
        results = [{"scientificName": "Homo sapiens", "synonyms": []}]
        assert bcfg._pick_taxonomy_match("Homo sapiens", results)["scientificName"] == "Homo sapiens"

    def test_taxonomy_is_vertebrate_true(self):
        payload = {
            "results": [
                {
                    "scientificName": "Homo sapiens",
                    "lineage": [{"scientificName": "Metazoa"}, {"taxonId": 7742}],
                }
            ]
        }
        with patch.object(bcfg, "uniprot_get", return_value=payload):
            assert bcfg.taxonomy_is_vertebrate("Homo sapiens") is True

    def test_taxonomy_is_vertebrate_invertebrate(self):
        payload = {
            "results": [
                {
                    "scientificName": "Drosophila melanogaster",
                    "lineage": [{"scientificName": "Metazoa"}, {"scientificName": "Ecdysozoa"}],
                }
            ]
        }
        with patch.object(bcfg, "uniprot_get", return_value=payload):
            assert bcfg.taxonomy_is_vertebrate("Drosophila melanogaster") is False


class TestDetectFamily:
    def test_detect_family_vertebrate(self):
        with patch.object(bcfg, "taxonomy_is_vertebrate", return_value=True):
            with patch("time.sleep"):
                family, method = bcfg.detect_family("Homo sapiens", True, ("connexin", "innexin"))
        assert family == "connexin"
        assert method == "taxonomy"

    def test_detect_family_invertebrate(self):
        with patch.object(bcfg, "taxonomy_is_vertebrate", return_value=False):
            with patch("time.sleep"):
                family, method = bcfg.detect_family("Drosophila melanogaster", True, ("connexin", "innexin"))
        assert family == "innexin"
        assert method == "taxonomy"

    def test_detect_family_uniprot_single_hit(self):
        with patch.object(bcfg, "taxonomy_is_vertebrate", return_value=None):
            with patch.object(bcfg, "search_uniprot_count", return_value=3):
                with patch("time.sleep"):
                    family, method = bcfg.detect_family("Unknown sp", True, ("innexin",))
        assert family == "innexin"
        assert method == "uniprot"

    def test_detect_family_default_fallback(self):
        with patch.object(bcfg, "taxonomy_is_vertebrate", return_value=None):
            with patch.object(bcfg, "search_uniprot_count", return_value=0):
                with patch("time.sleep"):
                    family, method = bcfg.detect_family("Unknown sp", True, ("connexin", "innexin"))
        assert family == "connexin"
        assert method == "default_fallback"


class TestConfigOutput:
    def test_build_protein_query(self):
        q = bcfg.build_protein_query("connexin", "Mus musculus", True)
        assert "protein_name:connexin" in q
        assert "reviewed:true" in q

    def test_write_config_csv(self, tmp_path: Path):
        out = tmp_path / "out.csv"
        rows = [{"family": "connexin", "organism": "Homo sapiens", "reviewed_only": "true"}]
        bcfg.write_config_csv(rows, out)
        with out.open(encoding="utf-8") as handle:
            parsed = list(csv.DictReader(handle))
        assert parsed[0]["organism"] == "Homo sapiens"

    def test_parse_args(self):
        args = bcfg.parse_args(["Homo sapiens", "-o", "out.csv"])
        assert args.species == ["Homo sapiens"]
        assert str(args.output) == "out.csv"
