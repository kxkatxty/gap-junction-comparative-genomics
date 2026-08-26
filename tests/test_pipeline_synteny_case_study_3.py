from __future__ import annotations

import json
from pathlib import Path

import pipeline.synteny_case_study_3_cross_clade as cs3


class TestSyntenyCaseStudy3:
    def test_build_clade_comparison_fly_only(self, tmp_path: Path, monkeypatch):
        synvoy = tmp_path / "SynVoy"
        fly_job = synvoy / "results" / "innexin_synvoy" / "dmel_inx2"
        (fly_job / "downloaded_genomes/easy_mode_genomes").mkdir(parents=True)
        (fly_job / "downloaded_genomes/easy_mode_genomes/species_mapping.tsv").write_text(
            "GCA_FLY.1\tDrosophila test\n",
            encoding="utf-8",
        )
        (fly_job / "synvoy_report.json").write_text(
            json.dumps(
                {
                    "synteny_results": {"genes_discovered": {"GCA_FLY.1": 5}},
                    "annotations": {
                        "per_genome": [
                            {
                                "genome": "GCA_FLY.1",
                                "goi_annotations": 1,
                                "role_counts": {"flanking": 4, "goi": 1},
                                "goi_class_counts": {"probable_goi": 1},
                                "goi_confidence_counts": {"LOW": 1},
                            }
                        ]
                    },
                    "regions": {"per_genome": [{"genome": "GCA_FLY.1", "total_regions": 1, "best_score": 0.2}]},
                }
            ),
            encoding="utf-8",
        )
        monkeypatch.setattr(cs3, "SYNVOY_DIR", synvoy)
        monkeypatch.setattr(cs3, "RESULTS_ROOT", fly_job.parent)
        import pipeline.synteny_case_study_2_ortholog as cs2

        monkeypatch.setattr(cs2, "SYNVOY_DIR", synvoy)
        monkeypatch.setattr(cs2, "RESULTS_ROOT", fly_job.parent)
        rows = cs3.build_clade_comparison()
        assert len(rows) == 2
        fly = next(r for r in rows if r["clade"] == "insect")
        cele = next(r for r in rows if r["clade"] == "nematode")
        assert fly["probable_ortholog_species"] == "1"
        assert cele["genomes_searched"] == "0"
