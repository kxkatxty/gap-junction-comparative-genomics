from __future__ import annotations

import extract_exon_structures as ex


class TestNormalization:
    def test_normalize_species_name(self):
        assert ex.normalize_species_name("Homo sapiens") == "Homo_sapiens"

    def test_normalize_gene_symbol(self):
        assert ex.normalize_gene_symbol("Gja-1") == "gja1"
        assert ex.normalize_gene_symbol("") == ""


class TestParseAttributes:
    def test_gff3_equals_format(self):
        attrs = ex.parse_attributes('ID=gene1;Name=GJA1;gene_biotype=protein_coding')
        assert attrs["ID"] == "gene1"
        assert attrs["Name"] == "GJA1"

    def test_gtf_quoted_format(self):
        attrs = ex.parse_attributes('gene_id "gene1"; transcript_id "tx1"')
        assert attrs["gene_id"] == "gene1"
        assert attrs["transcript_id"] == "tx1"


class TestReferenceParsing:
    def test_parse_reference_fasta(self, sample_fasta):
        ref = ex.parse_reference_fasta(sample_fasta)
        assert ref["gene_symbol"] == "GJA1"
        assert ref["uniprot_accession"] == "P17302"
        assert "Gap junction alpha-1" in ref["protein_name"]

    def test_candidate_gene_symbols_includes_aliases(self):
        ref = {
            "gene_symbol": "shakB",
            "uniprot_accession": "Q1DH70",
            "fasta_header": ">sp|Q1DH70|SHAK_DROME shaking-B OS=Drosophila GN=shakB",
            "protein_name": "Protein shakB",
        }
        symbols = ex.candidate_gene_symbols(ref)
        norms = {ex.normalize_gene_symbol(s) for s in symbols}
        assert "shakb" in norms
        assert "shakingb" in norms or "shaking-b" in symbols


class TestGeneMatching:
    def _build_indexes(self, mini_gff3):
        genes, transcripts, exons, gene_to_tx = ex.parse_annotation(mini_gff3)
        gene_index = ex.build_gene_index(genes)
        tx_index = ex.build_transcript_product_index(transcripts, genes)
        return genes, gene_index, tx_index, gene_to_tx

    def test_find_gja1_by_symbol(self, mini_gff3):
        genes, gene_index, tx_index, gene_to_tx = self._build_indexes(mini_gff3)
        reference = {
            "gene_symbol": "GJA1",
            "uniprot_accession": "P17302",
            "protein_name": "Gap junction alpha-1 protein",
            "fasta_header": "",
        }
        gene_id, method = ex.find_gene_for_reference(
            reference, genes, gene_index, tx_index, gene_to_tx
        )
        assert gene_id == "gene1"
        assert method.startswith("gene_symbol")

    def test_excludes_cnst_genes(self, mini_gff3):
        genes, _, _, _ = self._build_indexes(mini_gff3)
        cnst_gene = {
            "gene_symbol": "cnst",
            "gene_name": "connexin-associated protein",
        }
        assert ex.is_excluded_gene(cnst_gene)
        assert ex.is_excluded_gene(genes["gene2"]) is False

    def test_description_match_fallback(self, mini_gff3):
        genes, gene_index, tx_index, gene_to_tx = self._build_indexes(mini_gff3)
        reference = {
            "gene_symbol": "UNKNOWN",
            "uniprot_accession": "X12345",
            "protein_name": "Gap junction beta-2 protein",
            "fasta_header": "",
        }
        gene_id, method = ex.find_gene_for_reference(
            reference, genes, gene_index, tx_index, gene_to_tx
        )
        assert gene_id == "gene3"
        assert method == "description_match"


class TestTranscriptMetrics:
    def test_compute_transcript_metrics(self):
        exons = [(200, 400), (800, 1200), (3000, 4900)]
        metrics = ex.compute_transcript_metrics(exons)
        assert metrics["exon_count"] == 3
        assert metrics["total_exon_bases"] == 201 + 401 + 1901
        assert metrics["intron_count"] == 2
        assert metrics["transcript_span"] == 4900 - 200 + 1

    def test_empty_exons(self):
        metrics = ex.compute_transcript_metrics([])
        assert metrics["exon_count"] == 0
        assert metrics["transcript_span"] == 0


class TestGeneSymbolVariants:
    def test_variants_strip_prefix(self):
        variants = ex.gene_symbol_variants("GJA1")
        assert "gja1" in variants

    def test_product_search_keys(self):
        keys = ex.product_search_keys("Gap junction alpha-1 protein")
        assert "gap" in keys or "gapjunctionalpha1protein" in keys


class TestFamilyFolders:
    def test_family_to_ref_folder(self):
        assert ex.family_to_ref_folder("connexin").name == "connexins"
        assert ex.family_to_ref_folder("innexin").name == "innexins"


class TestUtilities:
    def test_safe_int(self):
        assert ex.safe_int("42") == 42
        assert ex.safe_int("3.0") == 3

    def test_load_species_csv(self, tmp_path):
        csv_path = tmp_path / "species.csv"
        csv_path.write_text(
            "family,organism,reviewed_only\nconnexin,Homo sapiens,true\n",
            encoding="utf-8",
        )
        rows = ex.load_species_csv(csv_path)
        assert rows[0]["family"] == "connexin"

    def test_read_fasta_header(self, sample_fasta):
        header = ex.read_fasta_header(sample_fasta)
        assert header.startswith(">")

    def test_description_match_score(self):
        reference = {"protein_name": "Gap junction alpha-1 protein"}
        gene = {"gene_name": "gap junction alpha-1 protein", "gene_symbol": "GJA1"}
        assert ex.description_match_score(reference, gene) >= 2

    def test_choose_best_gene_match_prefers_protein_coding(self):
        genes = {
            "g1": {"gene_symbol": "GJA1", "gene_biotype": "lncRNA"},
            "g2": {"gene_symbol": "GJA1", "gene_biotype": "protein_coding"},
        }
        best = ex.choose_best_gene_match(["g1", "g2"], genes, {"g2": ["tx1"]})
        assert best == "g2"

    def test_build_gene_index_excludes_cnst(self):
        genes = {
            "g1": {
                "gene_symbol": "cnst",
                "gene_name": "connexin-associated",
                "gene_synonym": "",
                "gene_biotype": "protein_coding",
            }
        }
        index = ex.build_gene_index(genes)
        assert index == {}

    def test_cxa_to_gja_mapping(self):
        ref = {
            "gene_symbol": "",
            "uniprot_accession": "P17302",
            "fasta_header": ">sp|P17302|CXA1_HUMAN",
            "protein_name": "",
        }
        symbols = ex.candidate_gene_symbols(ref)
        norms = {ex.normalize_gene_symbol(s) for s in symbols}
        assert "gja1" in norms

