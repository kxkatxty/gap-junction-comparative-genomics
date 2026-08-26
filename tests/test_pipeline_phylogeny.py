from __future__ import annotations

import csv
from pathlib import Path
from unittest.mock import MagicMock, patch

import pipeline.phylogeny_input_builder as phylo_input
import pipeline.phylogeny_runner as phylo_runner
import pipeline.sequence_alignment_runner as align_runner
import pipeline.tree_annotation_helper as tree_helper


class TestPhylogenyInputBuilder:
    def test_read_fasta_records(self, tmp_path: Path):
        fasta = tmp_path / "seq.fasta"
        fasta.write_text(">hdr1\nACGT\n>hdr2\nTGCA\n", encoding="utf-8")
        records = phylo_input.read_fasta_records(fasta)
        assert records == [("hdr1", "ACGT"), ("hdr2", "TGCA")]

    def test_build_records_references_only(self, pipeline_sandbox):
        records = phylo_input.build_records(
            family="innexin",
            include_references=True,
            accepted_only=True,
        )
        assert len(records) == 2
        labels = [r[0] for r in records]
        assert any("Q9V427" in label for label in labels)
        assert any("candidate_inx_001" in label for label in labels)

    def test_main_writes_fasta(self, pipeline_sandbox, tmp_path: Path):
        out = pipeline_sandbox["results"] / "phylo.fasta"
        assert phylo_input.main(["--family", "innexin", "--include-references", "--output", str(out)]) == 0
        text = out.read_text(encoding="utf-8")
        assert text.startswith(">")
        assert "MKTAYIAKQRQISFVKSHFSRQ" in text.replace("\n", "")

    def test_main_fails_without_sequences(self, tmp_path: Path, monkeypatch):
        empty_results = tmp_path / "results"
        empty_refs = tmp_path / "refs"
        empty_results.mkdir()
        empty_refs.mkdir()
        monkeypatch.setattr(phylo_input, "RESULTS_DIR", empty_results)
        monkeypatch.setattr(phylo_input, "REFERENCES_DIR", empty_refs)
        monkeypatch.setattr(phylo_input, "collect_candidates", lambda **kwargs: [])
        out = tmp_path / "empty.fasta"
        assert phylo_input.main(["--output", str(out)]) == 1


class TestTreeAnnotationHelper:
    def test_parse_newick_tips(self, tmp_path: Path):
        tree = tmp_path / "tree.nwk"
        tree.write_text("((A_species|x:0.1,B_species|y:0.2):0.3,C_species|z:0.4);", encoding="utf-8")
        tips = tree_helper.parse_newick_tips(tree)
        assert "A_species|x" in tips
        assert "C_species|z" in tips

    def test_build_labels_uses_master_summary(self, pipeline_sandbox, tmp_path: Path):
        master = pipeline_sandbox["metadata"] / "thesis_master_summary.csv"
        summary_rows = [
            {
                "family": "innexin",
                "organism": "Drosophila_melanogaster",
                "panel_role": "reference_innexin",
                "annotation_category": "annotated_family_member_found",
            }
        ]
        with master.open("w", encoding="utf-8", newline="") as handle:
            writer = csv.DictWriter(handle, fieldnames=summary_rows[0].keys())
            writer.writeheader()
            writer.writerows(summary_rows)

        tree = tmp_path / "tree.nwk"
        tree.write_text("(Drosophila_melanogaster|Inx2:0.1);", encoding="utf-8")
        rows = tree_helper.build_labels(tree)
        assert rows[0]["organism"] == "Drosophila_melanogaster"
        assert rows[0]["family"] == "innexin"
        assert "innexin" in rows[0]["itol_label"]


class TestSequenceAlignmentRunner:
    def test_main_missing_input(self, tmp_path: Path):
        assert align_runner.main([str(tmp_path / "missing.fasta")]) == 1

    @patch("pipeline.sequence_alignment_runner.shutil.which", return_value=None)
    def test_main_missing_mafft(self, mock_which, tmp_path: Path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">a\nACGT\n", encoding="utf-8")
        assert align_runner.main([str(fasta)]) == 1

    @patch("pipeline.sequence_alignment_runner.subprocess.run")
    @patch("pipeline.sequence_alignment_runner.shutil.which", return_value="/usr/bin/mafft")
    def test_main_runs_mafft(self, mock_which, mock_run, tmp_path: Path):
        fasta = tmp_path / "in.fasta"
        fasta.write_text(">a\nACGT\n", encoding="utf-8")
        mock_run.return_value = MagicMock(stdout=">a\nACGT\n")
        assert align_runner.main([str(fasta)]) == 0
        mock_run.assert_called_once()
        assert fasta.with_suffix(".aln.fasta").exists()


class TestPhylogenyRunner:
    def test_main_missing_alignment(self, tmp_path: Path):
        assert phylo_runner.main([str(tmp_path / "missing.aln.fasta")]) == 1

    @patch("pipeline.phylogeny_runner.shutil.which", return_value=None)
    def test_main_missing_iqtree(self, mock_which, tmp_path: Path):
        aln = tmp_path / "in.aln.fasta"
        aln.write_text(">a\nACGT\n", encoding="utf-8")
        assert phylo_runner.main([str(aln)]) == 1

    @patch("pipeline.phylogeny_runner.subprocess.run")
    def test_main_runs_iqtree(self, mock_run, tmp_path: Path, monkeypatch):
        import pipeline.phylogeny_runner as pr

        monkeypatch.setattr(pr, "DEFAULT_OUT_DIR", tmp_path)
        aln = tmp_path / "in.aln.fasta"
        aln.write_text(">a\nACGT\n", encoding="utf-8")

        with patch("pipeline.phylogeny_runner.shutil.which") as mock_which:
            mock_which.side_effect = lambda name: "/usr/bin/iqtree" if name == "iqtree" else None
            assert phylo_runner.main([str(aln), "--threads", "1"]) == 0

        mock_run.assert_called_once()
        cmd = mock_run.call_args[0][0]
        assert cmd[0] in {"iqtree", "iqtree2"}
        assert "-s" in cmd
