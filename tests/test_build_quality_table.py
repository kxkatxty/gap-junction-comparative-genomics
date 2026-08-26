from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import build_quality_table as bqt


class TestAccessionExtraction:
    def test_extract_from_filename(self, sample_fasta):
        header = sample_fasta.read_text(encoding="utf-8").splitlines()[0]
        assert bqt.extract_uniprot_accession(sample_fasta, header) == "P17302"

    def test_extract_from_header_sp_format(self, tmp_path: Path):
        path = tmp_path / "protein.fasta"
        path.write_text(">sp|P17302|CXA1_HUMAN Gap junction\nACDE\n", encoding="utf-8")
        header = bqt.read_first_fasta_header(path)
        assert bqt.extract_uniprot_accession(path, header) == "P17302"

    def test_extract_missing_returns_none(self, tmp_path: Path):
        path = tmp_path / "bad.fasta"
        path.write_text(">unknown header\nACDE\n", encoding="utf-8")
        header = bqt.read_first_fasta_header(path)
        assert bqt.extract_uniprot_accession(path, header) is None


class TestParseUniProt:
    def test_parse_uniprot_full_entry(self):
        entry = {
            "primaryAccession": "P17302",
            "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "organism": {"scientificName": "Homo sapiens"},
            "genes": [{"geneName": {"value": "GJA1"}}],
            "proteinDescription": {
                "recommendedName": {"fullName": {"value": "Gap junction alpha-1 protein"}}
            },
            "sequence": {"length": 382},
            "proteinExistence": "1: Evidence at protein level",
            "annotationScore": 5,
        }
        parsed = bqt.parse_uniprot(entry)
        assert parsed["accession"] == "P17302"
        assert parsed["gene_name"] == "GJA1"
        assert parsed["review_status"] == "reviewed"
        assert parsed["sequence_length"] == 382
        assert parsed["protein_existence"] == "1"

    def test_parse_uniprot_submission_name_fallback(self):
        entry = {
            "primaryAccession": "A0A000",
            "entryType": "UniProtKB TrEMBL",
            "organism": {"scientificName": "Test"},
            "genes": [],
            "proteinDescription": {
                "submissionNames": [{"fullName": {"value": "Submitted name"}}]
            },
            "sequence": {},
        }
        parsed = bqt.parse_uniprot(entry)
        assert parsed["protein_name"] == "Submitted name"
        assert parsed["review_status"] == "unreviewed"


class TestNcbiAssembly:
    def test_parse_ncbi_assembly_empty(self):
        parsed = bqt.parse_ncbi_assembly(None)
        assert parsed["assembly_accession"] == ""

    def test_parse_ncbi_assembly_refseq(self):
        summary = {
            "assemblyaccession": "GCF_000001405.40",
            "assemblystatus": "Chromosome",
            "assemblyname": "GRCh38.p14",
            "rsuid": "GCF_000001405.40",
            "meta": {"contig_n50": "100", "scaffold_n50": "200"},
        }
        parsed = bqt.parse_ncbi_assembly(summary)
        assert parsed["refseq_or_genbank"] == "RefSeq"
        assert parsed["assembly_level"] == "Chromosome"
        assert parsed["annotation_available"] == "yes"

    def test_ncbi_search_assembly_picks_best(self):
        search_payload = {"esearchresult": {"idlist": ["1", "2"]}}
        summary_payload = {
            "result": {
                "uids": ["1", "2"],
                "1": {"assemblyaccession": "GCA_1", "assemblystatus": "Contig", "rsuid": ""},
                "2": {
                    "assemblyaccession": "GCF_2",
                    "assemblystatus": "Chromosome",
                    "rsuid": "GCF_2",
                    "meta": {"scaffold_n50": "5000000"},
                },
            }
        }
        with patch.object(bqt, "http_get_json", side_effect=[search_payload, summary_payload]):
            best = bqt.ncbi_search_assembly("Homo sapiens")
        assert best["assemblyaccession"] == "GCF_2"


class TestQualityFlags:
    def test_make_quality_flag_high(self):
        assert bqt.make_quality_flag("reviewed", "Chromosome", "5") == "high"

    def test_make_quality_flag_medium(self):
        assert bqt.make_quality_flag("reviewed", "Scaffold", "2") == "medium"

    def test_make_quality_flag_review_protein_only(self):
        assert bqt.make_quality_flag("reviewed", "Contig", "1") == "review_protein_only"

    def test_make_quality_flag_low(self):
        assert bqt.make_quality_flag("unreviewed", "Contig", "1") == "low_or_manual_check"


class TestHtmlReport:
    def test_write_html_report(self, tmp_path: Path):
        rows = [
            {
                "family": "connexin",
                "accession": "P17302",
                "organism": "Homo sapiens",
                "quality_flag": "high",
                "protein_name": "Gap junction alpha-1 protein",
                "gene_name": "GJA1",
                "review_status": "reviewed",
                "assembly_level": "Chromosome",
                "annotation_score": "5",
            }
        ]
        out = tmp_path / "report.html"
        bqt.write_html_report(rows, out)
        html = out.read_text(encoding="utf-8")
        assert "P17302" in html
        assert "Homo sapiens" in html
        assert "<html" in html.lower()
