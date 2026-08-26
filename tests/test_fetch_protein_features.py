from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import fetch_protein_features as fpf


class TestWhitelist:
    def test_is_enabled(self):
        assert fpf._is_enabled("true")
        assert fpf._is_enabled("1")
        assert fpf._is_enabled("yes")
        assert not fpf._is_enabled("false")
        assert not fpf._is_enabled(None)

    def test_match_whitelist_entries(self, whitelist_csv, mapping_csv, tmp_path):
        mapping = fpf.load_mapping(mapping_csv)
        mapping["fasta_path"] = mapping["fasta_path"].apply(
            lambda p: str(tmp_path / Path(p).name)
        )
        whitelist = fpf.load_whitelist(whitelist_csv)
        hits = fpf.match_whitelist_entries(whitelist, mapping)
        assert not hits.empty
        assert set(hits["plot_group"]) == {"GJA1_human", "GJA1"}
        human = hits[hits["plot_group"] == "GJA1_human"]
        assert len(human) == 1
        assert human.iloc[0]["uniprot_accession"] == "P17302"

    def test_whitelist_missing_columns(self, tmp_path):
        bad = tmp_path / "bad.csv"
        bad.write_text("plot_group,family\nGJA1,connexin\n", encoding="utf-8")
        with __import__("pytest").raises(ValueError, match="missing columns"):
            fpf.load_whitelist(bad)


class TestPfamParsing:
    def test_fetch_pfam_features(self):
        payload = {
            "results": [
                {
                    "metadata": {"accession": "PF00822", "name": "Connexin"},
                    "proteins": [
                        {
                            "entry_protein_locations": [
                                {"fragments": [{"start": 10, "end": 200}]}
                            ]
                        }
                    ],
                }
            ]
        }
        with patch.object(fpf, "http_get_json", return_value=payload):
            features = fpf.fetch_pfam_features("P17302")
        assert len(features) == 1
        assert features[0]["track"] == "pfam"
        assert features[0]["aa_start"] == 10
        assert features[0]["feature_id"] == "PF00822"

    def test_fetch_pfam_http_error(self):
        import urllib.error

        with patch.object(
            fpf,
            "http_get_json",
            side_effect=urllib.error.HTTPError(None, 404, "", {}, None),
        ):
            assert fpf.fetch_pfam_features("MISSING") == []


class TestUniprotParsing:
    def test_fetch_uniprot_coiled_and_disorder(self):
        payload = {
            "features": [
                {
                    "type": "Coiled coil",
                    "description": "",
                    "location": {"start": {"value": 5}, "end": {"value": 30}},
                },
                {
                    "type": "Region",
                    "description": "Intrinsically disordered region",
                    "location": {"start": {"value": 100}, "end": {"value": 150}},
                },
                {
                    "type": "Signal peptide",
                    "description": "",
                    "location": {"start": {"value": 1}, "end": {"value": 20}},
                },
            ]
        }
        with patch.object(fpf, "http_get_json", return_value=payload):
            features = fpf.fetch_uniprot_features("P17302")
        tracks = {f["track"] for f in features}
        assert "coiled_coil" in tracks
        assert "disorder" in tracks
        assert "signal" not in tracks

    def test_fetch_features_for_accession_combines(self):
        with patch.object(fpf, "fetch_pfam_features", return_value=[{"track": "pfam"}]):
            with patch.object(
                fpf, "fetch_uniprot_features", return_value=[{"track": "disorder"}]
            ):
                with patch("time.sleep"):
                    combined = fpf.fetch_features_for_accession("P17302")
        assert len(combined) == 2

    def test_load_mapping_missing_file(self, tmp_path):
        with __import__("pytest").raises(FileNotFoundError):
            fpf.load_mapping(tmp_path / "missing.csv")

