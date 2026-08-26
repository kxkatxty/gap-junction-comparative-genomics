from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import pipeline.case_study_synteny_runner as synteny
import pipeline.figure_export_pipeline as figures


class TestFigureExportPipeline:
    def test_collect_sources_finds_png_and_html(self, pipeline_sandbox):
        sources = figures.collect_sources()
        names = {name for _, name in sources}
        assert any(name.startswith("exon_") for name in names)
        assert any(name.startswith("synvoy_") for name in names)

    def test_main_exports_files(self, pipeline_sandbox):
        out = pipeline_sandbox["results"] / "figures" / "thesis"
        assert figures.main(["--output-dir", str(out), "--clean"]) == 0
        assert (out / "manifest.txt").exists()
        exported = list(out.glob("*"))
        assert len(exported) >= 2


class TestCaseStudySyntenyRunner:
    def test_load_queue(self, pipeline_sandbox):
        queue = pipeline_sandbox["metadata"] / "queue.tsv"
        queue.write_text(
            "gene\taccession\toutdir_name\n"
            "Inx2\tQ9V427\tdmel_inx2\n",
            encoding="utf-8",
        )
        rows = synteny.load_queue(queue)
        assert rows[0]["gene"] == "Inx2"

    def test_job_done_false_without_plot(self, pipeline_sandbox, monkeypatch):
        monkeypatch.setattr(synteny, "SYNVOY_DIR", pipeline_sandbox["synvoy"])
        assert synteny.job_done("missing_job") is False

    def test_job_done_true_with_plot(self, pipeline_sandbox, monkeypatch):
        monkeypatch.setattr(synteny, "SYNVOY_DIR", pipeline_sandbox["synvoy"])
        assert synteny.job_done("dmel_ogre") is True

    @patch("pipeline.case_study_synteny_runner.subprocess.run")
    def test_run_job_skips_completed(self, mock_run, pipeline_sandbox, monkeypatch):
        monkeypatch.setattr(synteny, "SYNVOY_DIR", pipeline_sandbox["synvoy"])
        row = {"gene": "ogre", "accession": "P27716", "outdir_name": "dmel_ogre"}
        assert synteny.run_job(row, max_genomes=3) == 0
        mock_run.assert_not_called()

    @patch("pipeline.case_study_synteny_runner.subprocess.run")
    def test_run_job_launches_nextflow(self, mock_run, pipeline_sandbox, monkeypatch):
        monkeypatch.setattr(synteny, "SYNVOY_DIR", pipeline_sandbox["synvoy"])
        monkeypatch.setattr(synteny, "LOG_DIR", pipeline_sandbox["results"] / "logs")
        mock_run.return_value = MagicMock(returncode=0)
        row = {"gene": "NewGene", "accession": "P00000", "outdir_name": "new_job"}
        assert synteny.run_job(row, max_genomes=2) == 0
        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] == "nextflow"
        assert "--query_id" in cmd
        assert "P00000" in cmd

    def test_main_no_matching_jobs(self, pipeline_sandbox, monkeypatch):
        queue = pipeline_sandbox["metadata"] / "queue.tsv"
        queue.write_text("gene\taccession\toutdir_name\nInx2\tQ9V427\tdmel_inx2\n", encoding="utf-8")
        monkeypatch.setattr(synteny, "DEFAULT_QUEUE", queue)
        assert synteny.main(["--gene", "MissingGene"]) == 1
