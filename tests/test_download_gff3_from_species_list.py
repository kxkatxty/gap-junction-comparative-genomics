from __future__ import annotations

import gzip
from pathlib import Path
from unittest.mock import MagicMock, patch

import download_gff3_from_species_list as dgff


class TestSanitizeAndCsv:
    def test_sanitize_filename(self):
        assert dgff.sanitize_filename("Homo sapiens (GRCh38)") == "Homo_sapiens_GRCh38"

    def test_read_species_csv(self, tmp_path: Path):
        csv_path = tmp_path / "species.csv"
        csv_path.write_text(
            "family,organism,reviewed_only\nconnexin,Homo sapiens,true\n",
            encoding="utf-8",
        )
        rows = dgff.read_species_csv(csv_path)
        assert rows[0]["organism"] == "Homo sapiens"

    def test_read_species_csv_missing_columns(self, tmp_path: Path):
        bad = tmp_path / "bad.csv"
        bad.write_text("family,organism\nconnexin,Homo sapiens\n", encoding="utf-8")
        with __import__("pytest").raises(ValueError, match="Missing required columns"):
            dgff.read_species_csv(bad)


class TestAssemblyRanking:
    def test_assembly_search_terms(self):
        terms = dgff.assembly_search_terms("Homo sapiens")
        assert any("latest_refseq" in t for t in terms)
        assert all("Homo sapiens" in t for t in terms)

    def test_rank_assemblies_prefers_refseq_chromosome(self):
        records = [
            {
                "AssemblyAccession": "GCA_000001",
                "FtpPath_RefSeq": "",
                "RefSeq_category": "",
                "AssemblyStatus": "Scaffold",
            },
            {
                "AssemblyAccession": "GCF_000002",
                "FtpPath_RefSeq": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/002",
                "RefSeq_category": "reference genome",
                "AssemblyStatus": "Chromosome",
            },
        ]
        ranked = dgff.rank_assemblies(records)
        assert ranked[0]["AssemblyAccession"] == "GCF_000002"

    def test_rank_assemblies_deduplicates_by_accession(self):
        records = [
            {"AssemblyAccession": "GCF_1", "FtpPath_RefSeq": "", "AssemblyStatus": "Contig"},
            {"AssemblyAccession": "GCF_1", "FtpPath_RefSeq": "ftp://x", "AssemblyStatus": "Chromosome"},
        ]
        ranked = dgff.rank_assemblies(records)
        assert len(ranked) == 1
        assert "Chromosome" in ranked[0]["AssemblyStatus"]

    def test_gff_candidate_urls(self):
        rec = {"FtpPath_RefSeq": "ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCF/000/001"}
        urls = dgff.gff_candidate_urls(rec)
        assert urls[0][2] == "RefSeq"
        assert urls[0][0].startswith("https://")
        assert "genomic.gff.gz" in urls[0][0]


class TestEntrezParsing:
    def test_entrez_esearch_assembly(self):
        xml = """<?xml version="1.0"?>
        <eSearchResult><IdList><Id>123</Id><Id>456</Id></IdList></eSearchResult>"""
        with patch.object(dgff, "http_get_text", return_value=xml):
            ids = dgff.entrez_esearch_assembly("Homo sapiens", "term")
        assert ids == ["123", "456"]

    def test_entrez_esummary_assembly(self):
        xml = """<?xml version="1.0"?>
        <eSummaryResult>
          <DocumentSummary uid="1"><AssemblyAccession>GCF_000</AssemblyAccession></DocumentSummary>
        </eSummaryResult>"""
        with patch.object(dgff, "http_get_text", return_value=xml):
            with patch("time.sleep"):
                records = dgff.entrez_esummary_assembly(["1"])
        assert records[0]["AssemblyAccession"] == "GCF_000"


class TestFileHelpers:
    def test_existing_gff_files(self, tmp_path: Path):
        (tmp_path / "asm_genomic.gff3").write_text("##gff\n", encoding="utf-8")
        found = dgff.existing_gff_files(tmp_path)
        assert len(found) == 1

    def test_maybe_decompress_gzip(self, tmp_path: Path):
        gz_path = tmp_path / "test.gff.gz"
        with gzip.open(gz_path, "wb") as handle:
            handle.write(b"##gff-version 3\n")
        out = dgff.maybe_decompress_gzip(gz_path, keep_gz=True)
        assert out.exists()
        assert out.read_text(encoding="utf-8").startswith("##gff")

    def test_download_file(self, tmp_path: Path):
        out = tmp_path / "out.gff"
        payload = b"gff content" + (b"x" * 1200)
        with patch.object(dgff, "urlopen") as mock_open, patch.object(dgff, "_validate_downloaded_file"):
            mock_resp = MagicMock()
            remaining = {"data": payload}

            def read(size=-1):
                if not remaining["data"]:
                    return b""
                if size is None or size < 0:
                    chunk = remaining["data"]
                else:
                    chunk = remaining["data"][:size]
                remaining["data"] = remaining["data"][len(chunk) :]
                return chunk

            mock_resp.read = read
            mock_resp.__enter__.return_value = mock_resp
            mock_open.return_value = mock_resp
            dgff.download_file("https://example.com/a.gff", out)
        assert out.read_bytes() == payload


class TestDownloadGffForOrganism:
    def test_skips_when_gff_exists(self, tmp_path: Path):
        organism_dir = tmp_path / "Homo_sapiens"
        organism_dir.mkdir()
        (organism_dir / "existing.gff3").write_text("##gff\n", encoding="utf-8")
        result = dgff.download_gff_for_organism("Homo sapiens", organism_dir, decompress=False)
        assert result["status"] == "skipped_existing"

    def test_raises_when_no_candidates(self, tmp_path: Path):
        organism_dir = tmp_path / "Unknown_sp"
        organism_dir.mkdir()
        with patch.object(dgff, "fetch_assembly_candidates", return_value=[]):
            with __import__("pytest").raises(ValueError, match="No assembly candidates"):
                dgff.download_gff_for_organism("Unknown sp", organism_dir, decompress=False)


class TestHttpHelpers:
    def test_url_exists_true(self):
        with patch.object(dgff, "urlopen") as mock_open:
            mock_resp = MagicMock()
            mock_resp.status = 200
            mock_resp.__enter__.return_value = mock_resp
            mock_open.return_value = mock_resp
            assert dgff.url_exists("https://example.com/file.gff")

    def test_url_exists_404(self):
        import urllib.error

        with patch.object(dgff, "urlopen") as mock_open:
            mock_open.side_effect = urllib.error.HTTPError(
                "url", 404, "", {}, None
            )
            assert dgff.url_exists("https://example.com/missing.gff") is False

    def test_parse_args_defaults(self):
        args = dgff.parse_args(["species_config.csv"])
        assert args.csv == "species_config.csv"
        assert args.no_decompress is False
