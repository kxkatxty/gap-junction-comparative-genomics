from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import fetch_biocentral_predictions as fbp


class TestSequenceUtils:
    def test_read_fasta_sequence(self, sample_fasta):
        assert fbp.read_fasta_sequence(sample_fasta) == "MKTAYIAKQRQISFVKSHFSRQ"

    def test_sequence_hash_stable(self):
        assert fbp.sequence_hash("ACDE") == fbp.sequence_hash("ACDE")
        assert len(fbp.sequence_hash("ACDE")) == 16


class TestIntervalMerging:
    def test_merge_class_intervals_helix(self):
        values = "LLLHHHHLLL"
        class_map = {
            "H": ("helix", "Helix", "#FF0000"),
            "L": ("coil", "Coil", "#CCCCCC"),
        }
        features = fbp.merge_class_intervals(
            values, class_map, "P1", "secondary_structure", "test"
        )
        helix = [f for f in features if f["feature_id"] == "helix"]
        assert len(helix) == 1
        assert helix[0]["aa_start"] == 4
        assert helix[0]["aa_end"] == 7

    def test_merge_class_intervals_skips_short_regions(self):
        values = "HH"
        class_map = {"H": ("helix", "Helix", "#FF0000")}
        features = fbp.merge_class_intervals(
            values, class_map, "P1", "secondary_structure", "test"
        )
        assert features == []

    def test_merge_threshold_intervals_disorder(self):
        scores = [10.0, 9.0, 5.0, 4.0, 4.0, 11.0]
        features = fbp.merge_threshold_intervals(
            scores, "P1", threshold=8.0, below_is_positive=True
        )
        assert len(features) == 1
        assert features[0]["aa_start"] == 3
        assert features[0]["aa_end"] == 5
        assert features[0]["source"] == "biocentral_seth"


class TestPredictionParsing:
    def test_parse_secondary_structure_excludes_coil(self):
        values = "LLLHHHHHLLL"
        features = fbp.parse_secondary_structure("P1", values)
        assert all(f["feature_id"] != "coil" for f in features)
        assert any(f["feature_id"] == "helix" for f in features)

    def test_parse_transmembrane(self):
        values = "iii" + "H" * 5 + "o" * 3
        features = fbp.parse_transmembrane("P1", values)
        assert len(features) == 1
        assert features[0]["feature_id"] == "tm_helix"
        assert features[0]["aa_start"] == 4

    def test_parse_seth_disorder_from_list(self):
        scores = [10, 10, 5, 4, 3, 10]
        features = fbp.parse_seth_disorder("P1", scores)
        assert len(features) == 1
        assert features[0]["track"] == "disorder"

    def test_predictions_to_features_full_payload(self):
        predictions = [
            {
                "model_name": "ProtT5SecondaryStructure",
                "prediction_name": "d3_Yhat",
                "value": "LLLHHHHLLL",
            },
            {
                "model_name": "TMbed",
                "prediction_name": "trans_membrane",
                "value": "iiiHHHHHooo",
            },
            {
                "model_name": "Seth",
                "prediction_name": "disorder_chezod",
                "value": [10, 10, 5, 4, 3, 10],
            },
        ]
        features = fbp.predictions_to_features("P17302", predictions)
        tracks = {f["track"] for f in features}
        assert "secondary_structure" in tracks
        assert "transmembrane" in tracks
        assert "disorder" in tracks


class TestApiWrappers:
    def test_submit_prediction_batch(self):
        with patch.object(fbp, "http_json", return_value={"task_id": "task-123"}):
            task_id = fbp.submit_prediction_batch({"P1": "ACDE"})
        assert task_id == "task-123"

    def test_poll_task_finished(self):
        dto = {"status": "FINISHED", "predictions": {}}
        with patch("time.sleep"):
            with patch.object(
                fbp,
                "http_json",
                return_value={"dtos": [dto]},
            ):
                result = fbp.poll_task("task-123")
        assert result["status"] == "FINISHED"

    def test_poll_task_retries_on_403(self):
        import urllib.error

        dto = {"status": "FINISHED", "predictions": {}}
        responses = [
            urllib.error.HTTPError("url", 403, "", {}, None),
            {"dtos": [dto]},
        ]

        def side_effect(*_args, **_kwargs):
            item = responses.pop(0)
            if isinstance(item, Exception):
                raise item
            return item

        with patch("time.sleep"):
            with patch.object(fbp, "http_json", side_effect=side_effect):
                result = fbp.poll_task("task-123")
        assert result["status"] == "FINISHED"


class TestCacheAndTargets:
    def test_save_cache_writes_json(self, tmp_path, monkeypatch):
        monkeypatch.setattr(fbp, "CACHE_DIR", tmp_path)
        fbp.save_cache("P17302", "task-1", {"status": "FINISHED"}, "abc123")
        cache_files = list(tmp_path.glob("P17302_*.json"))
        assert len(cache_files) == 1
        assert "task-1" in cache_files[0].read_text(encoding="utf-8")

    def test_load_existing_predictions_empty(self, tmp_path, monkeypatch):
        monkeypatch.setattr(fbp, "BIOCENTRAL_CSV", tmp_path / "missing.csv")
        df, cached = fbp.load_existing_predictions()
        assert df.empty
        assert cached == set()

    def test_load_existing_predictions_with_hash(self, tmp_path, monkeypatch):
        csv_path = tmp_path / "pred.csv"
        csv_path.write_text(
            "uniprot_accession,track,feature_id,aa_start,aa_end,sequence_hash\n"
            "P17302,disorder,disorder,1,10,hash1\n",
            encoding="utf-8",
        )
        monkeypatch.setattr(fbp, "BIOCENTRAL_CSV", csv_path)
        df, cached = fbp.load_existing_predictions()
        assert "hash1" in cached
        assert len(df) == 1

    def test_parse_seth_disorder_from_string(self):
        features = fbp.parse_seth_disorder("P1", "10,9,5,4,4,11")
        assert len(features) == 1

    def test_merge_threshold_intervals_above_threshold(self):
        scores = [1.0, 2.0, 10.0, 10.0, 10.0, 11.0]
        features = fbp.merge_threshold_intervals(
            scores, "P1", threshold=5.0, below_is_positive=False
        )
        assert len(features) == 1
        assert features[0]["aa_start"] == 3

