from __future__ import annotations

import json
from pathlib import Path

import pipeline.synteny_case_study_2_ortholog as cs2


class TestSyntenyCaseStudy2:
    def test_classify_support(self):
        assert cs2.classify_support({"confident_goi": "1"}) == "confident_ortholog"
        assert cs2.classify_support({"probable_goi": "1"}) == "probable_ortholog"
        assert cs2.classify_support({"goi_annotations": "2", "ambiguous_goi": "2"}) == "weak_goi_signal"
        assert cs2.classify_support({"flanking_annotations": "5"}) == "flanking_synteny_only"
        assert cs2.classify_support({}) == "no_support"

    def test_build_ortholog_rows(self, tmp_path: Path, monkeypatch):
        synvoy = tmp_path / "SynVoy"
        job = synvoy / "results" / "innexin_synvoy" / "dmel_inx2"
        (job / "downloaded_genomes/easy_mode_genomes").mkdir(parents=True)
        (job / "downloaded_genomes/easy_mode_genomes/species_mapping.tsv").write_text(
            "GCA_TEST.1\tTest species\n",
            encoding="utf-8",
        )
        (job / "synvoy_report.json").write_text(
            json.dumps(
                {
                    "synteny_results": {"genes_discovered": {"GCA_TEST.1": 10}},
                    "annotations": {
                        "per_genome": [
                            {
                                "genome": "GCA_TEST.1",
                                "total_annotations": 10,
                                "role_counts": {"flanking": 8, "goi": 2},
                                "goi_annotations": 2,
                                "goi_confidence_counts": {"LOW": 2},
                                "goi_class_counts": {"probable_goi": 1},
                            }
                        ]
                    },
                    "regions": {
                        "per_genome": [
                            {"genome": "GCA_TEST.1", "total_regions": 2, "goi_anchor_regions": 1, "best_score": 0.5}
                        ]
                    },
                }
            ),
            encoding="utf-8",
        )
        monkeypatch.setattr(cs2, "SYNVOY_DIR", synvoy)
        monkeypatch.setattr(cs2, "RESULTS_ROOT", job.parent)
        rows = cs2.build_ortholog_rows(cs2.CASE_GENES[0])
        assert rows[0]["target_species"] == "Test species"
        assert rows[0]["ortholog_support"] == "probable_ortholog"
