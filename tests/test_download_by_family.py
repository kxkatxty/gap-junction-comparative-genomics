from __future__ import annotations

from pathlib import Path
from unittest.mock import MagicMock, patch

import download_by_family as dlf


class TestSanitizeAndFolders:
    def test_sanitize_filename(self):
        assert dlf.sanitize_filename("Gap junction α-1") == "Gap_junction__-1"

    def test_family_to_folder_connexin(self):
        assert dlf.family_to_folder("connexin").name == "connexins"

    def test_family_to_folder_innexin(self):
        assert dlf.family_to_folder("innexins").name == "innexins"

    def test_family_to_folder_pannexin(self):
        assert dlf.family_to_folder("pannexin").name == "pannexins"

    def test_family_to_folder_unsupported(self):
        with __import__("pytest").raises(ValueError, match="Unsupported family"):
            dlf.family_to_folder("unknown")


class TestUniProtQuery:
    def test_build_query_reviewed(self):
        q = dlf.build_query("connexin", "Homo sapiens", reviewed_only=True)
        assert "protein_name:connexin" in q
        assert 'organism_name:"Homo sapiens"' in q
        assert "reviewed:true" in q

    def test_build_query_unreviewed(self):
        q = dlf.build_query("innexin", "Drosophila melanogaster", reviewed_only=False)
        assert "reviewed:true" not in q

    def test_uniprot_organism_names_deduplicates(self):
        with patch.dict(dlf.UNIPROT_ORGANISM_ALIASES, {"Homo sapiens": ["Homo sapiens", "human"]}):
            names = dlf.uniprot_organism_names("Homo sapiens")
        assert names[0] == "Homo sapiens"
        assert len(names) == len(set(n.casefold() for n in names))


class TestEntryParsing:
    def test_extract_entry_info(self):
        entry = {
            "primaryAccession": "P17302",
            "entryType": "UniProtKB reviewed (Swiss-Prot)",
            "organism": {"scientificName": "Homo sapiens"},
            "genes": [{"geneName": {"value": "GJA1"}}],
            "proteinDescription": {
                "recommendedName": {"fullName": {"value": "Gap junction alpha-1 protein"}}
            },
        }
        info = dlf.extract_entry_info(entry)
        assert info["accession"] == "P17302"
        assert info["gene_name"] == "GJA1"
        assert info["review_status"] == "reviewed"

    def test_extract_entry_info_unreviewed(self):
        entry = {
            "primaryAccession": "A0A000",
            "entryType": "UniProtKB unreviewed (TrEMBL)",
            "organism": {"scientificName": "Test sp"},
            "genes": [],
            "proteinDescription": {},
        }
        info = dlf.extract_entry_info(entry)
        assert info["review_status"] == "unreviewed"


class TestFastaIO:
    def test_save_fasta(self, tmp_path: Path):
        folder = tmp_path / "connexins"
        fasta = ">sp|P17302|CXA1_HUMAN test\nACDE\n"
        out = dlf.save_fasta(folder, "Homo sapiens", "GJA1", "P17302", fasta)
        assert out.exists()
        assert out.name == "GJA1__P17302.fasta"
        assert "P17302" in out.read_text(encoding="utf-8")

    def test_fetch_fasta_valid(self):
        with patch.object(dlf, "urlopen") as mock_open:
            mock_resp = MagicMock()
            mock_resp.read.return_value = b">sp|P17302|CXA1\nACDE\n"
            mock_resp.__enter__.return_value = mock_resp
            mock_open.return_value = mock_resp
            text = dlf.fetch_fasta("P17302")
        assert text.startswith(">")

    def test_fetch_fasta_invalid(self):
        with patch.object(dlf, "urlopen") as mock_open:
            mock_resp = MagicMock()
            mock_resp.read.return_value = b"not fasta"
            mock_resp.__enter__.return_value = mock_resp
            mock_open.return_value = mock_resp
            with __import__("pytest").raises(ValueError, match="FASTA"):
                dlf.fetch_fasta("BAD")

    def test_search_uniprot(self):
        payload = {"results": [{"primaryAccession": "P17302"}]}
        with patch.object(dlf, "urlopen") as mock_open:
            mock_resp = MagicMock()
            mock_resp.read.return_value = __import__("json").dumps(payload).encode()
            mock_resp.__enter__.return_value = mock_resp
            mock_open.return_value = mock_resp
            result = dlf.search_uniprot("connexin AND reviewed:true")
        assert result["results"][0]["primaryAccession"] == "P17302"
