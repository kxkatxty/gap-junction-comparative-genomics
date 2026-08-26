from __future__ import annotations

import csv
import json
from pathlib import Path

import pipeline.synteny_comparison_builder as scb


class TestSyntenyComparisonBuilder:
    def test_parse_home_flanking(self, tmp_path: Path):
        bed = tmp_path / "block.bed"
        bed.write_text(
            "chr1\t100\t200\tgene-A\t.\t+\n"
            "chr1\t300\t400\tGOI_chr1_350\t.\t-\n"
            "chr1\t500\t600\tgene-B\t.\t+\n",
            encoding="utf-8",
        )
        rows = scb.parse_home_flanking(bed)
        assert len(rows) == 3
        assert rows[1]["is_goi"] == "yes"
        assert rows[0]["is_goi"] == "no"

    def test_build_job_summary_from_report(self, tmp_path: Path, monkeypatch):
        synvoy = tmp_path / "SynVoy"
        results = synvoy / "results" / "innexin_synvoy" / "dmel_test"
        results.mkdir(parents=True)
        (results / "synteny_block_locus_1_synteny_plot.html").write_text("<html></html>", encoding="utf-8")
        (results / "synvoy_report.json").write_text(
            json.dumps(
                {
                    "qc_summary": {"total_genomes": 3},
                    "summary": {
                        "genomes_with_hits": 2,
                        "total_hits": 5,
                        "genomes_with_annotations": 2,
                        "total_goi_annotations": 4,
                        "low_confidence_regions": 1,
                    },
                    "annotations": {"per_genome": [{"goi_class_counts": {"confident_goi": 1}}]},
                }
            ),
            encoding="utf-8",
        )
        queue = tmp_path / "queue.tsv"
        queue.write_text(
            "accession\tgene\tsource_species\toutdir_name\n"
            "P00000\tInx2\tDrosophila_melanogaster\tdmel_test\n",
            encoding="utf-8",
        )
        monkeypatch.setattr(scb, "SYNVOY_DIR", synvoy)
        monkeypatch.setattr(scb, "RESULTS_ROOT", results.parent)

        rows = scb.build_job_summary(queue)
        assert len(rows) == 1
        assert rows[0]["status"] == "complete"
        assert rows[0]["gene"] == "Inx2"
        assert rows[0]["genomes_searched"] == "3"

    def test_main_writes_tables(self, tmp_path: Path, monkeypatch):
        synvoy = tmp_path / "SynVoy"
        results = synvoy / "results" / "innexin_synvoy" / "dmel_test"
        plot_inputs = results / "plot_inputs_synteny_block_locus_1"
        plot_inputs.mkdir(parents=True)
        (results / "synteny_block_locus_1_synteny_plot.html").write_text("<html></html>", encoding="utf-8")
        (plot_inputs / "synteny_block_locus_1.bed").write_text(
            "chr1\t100\t200\tflank1\t.\t+\nchr1\t300\t400\tGOI_chr1_350\t.\t-\n",
            encoding="utf-8",
        )
        (results / "synvoy_report.json").write_text(
            json.dumps(
                {
                    "qc_summary": {"total_genomes": 1},
                    "summary": {},
                    "synteny_results": {"genes_discovered": {"GCA_1.1": 10}, "synteny_hits_count": {}},
                    "annotations": {
                        "per_genome": [
                            {
                                "genome": "GCA_1.1",
                                "total_annotations": 10,
                                "role_counts": {"flanking": 8, "goi": 2},
                                "goi_annotations": 2,
                                "goi_class_counts": {"confident_goi": 1},
                                "goi_confidence_counts": {"HIGH": 1},
                            }
                        ]
                    },
                }
            ),
            encoding="utf-8",
        )
        queue = tmp_path / "queue.tsv"
        queue.write_text(
            "accession\tgene\tsource_species\toutdir_name\n"
            "P00000\tInx2\tDrosophila_melanogaster\tdmel_test\n",
            encoding="utf-8",
        )
        out_dir = tmp_path / "metadata"
        monkeypatch.setattr(scb, "SYNVOY_DIR", synvoy)
        monkeypatch.setattr(scb, "RESULTS_ROOT", results.parent)
        assert scb.main(["--queue", str(queue), "--output-dir", str(out_dir)]) == 0
        assert (out_dir / "synteny_job_summary.csv").exists()
        assert (out_dir / "synteny_genome_comparison.csv").exists()
        assert (out_dir / "synteny_home_flanking_genes.csv").exists()
        genome_rows = list(csv.DictReader((out_dir / "synteny_genome_comparison.csv").open(encoding="utf-8")))
        assert genome_rows[0]["confident_goi"] == "1"
