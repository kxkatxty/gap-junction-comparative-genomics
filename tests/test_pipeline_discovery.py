from __future__ import annotations

import csv

import pipeline.candidate_validator as validator
import pipeline.family_result_merger as merger


class TestFamilyResultMerger:
    def test_collect_candidates_all(self, pipeline_sandbox):
        rows = merger.collect_candidates(accepted_only=False)
        assert len(rows) == 2
        assert {r["rank_category"] for r in rows} == {
            "high_confidence_innexin_candidate",
            "rejected_false_positive",
        }

    def test_collect_candidates_accepted_only(self, pipeline_sandbox):
        rows = merger.collect_candidates(accepted_only=True)
        assert len(rows) == 1
        assert rows[0]["candidate_id"] == "candidate_inx_001"

    def test_main_writes_master_csv(self, pipeline_sandbox):
        out = pipeline_sandbox["metadata"] / "candidates.csv"
        assert merger.main(["--output", str(out)]) == 0
        rows = list(csv.DictReader(out.open(encoding="utf-8")))
        assert len(rows) == 2


class TestCandidateValidator:
    def test_classify_confidence(self):
        assert validator.classify_confidence("rejected_false_positive") == "reject"
        assert validator.classify_confidence("high_confidence_innexin_candidate") == "high"
        assert validator.classify_confidence("weak_manual_review") == "weak"
        assert validator.classify_confidence("possible_connexin_candidate") == "possible"

    def test_build_summary_counts(self, pipeline_sandbox):
        rows = validator.build_summary()
        assert len(rows) == 1
        row = rows[0]
        assert row["organism"] == "Test species"
        assert row["candidate_total"] == "2"
        assert row["accepted_total"] == "1"
        assert row["high_confidence"] == "1"
        assert row["rejected"] == "1"
        assert row["best_rank_category"] == "high_confidence_innexin_candidate"

    def test_main_writes_summary(self, pipeline_sandbox):
        out = pipeline_sandbox["metadata"] / "validation.csv"
        assert validator.main(["--output", str(out)]) == 0
        rows = list(csv.DictReader(out.open(encoding="utf-8")))
        assert rows[0]["accepted_total"] == "1"
