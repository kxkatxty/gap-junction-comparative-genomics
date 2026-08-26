from __future__ import annotations

from unittest.mock import patch

import classify_family_annotation as cfa


class TestClassificationLogic:
    def test_exact_family_word_is_annotated(self):
        hit = cfa.SearchHit(
            database="uniprot",
            level="A",
            query="q",
            protein_name="Gap junction alpha-1 protein (connexin)",
            gene_name="GJA1",
        )
        result = cfa.classify_hit("connexin", hit)
        assert result is not None
        assert result.final_category == cfa.CATEGORY_ANNOTATED
        assert result.confidence == cfa.CONFIDENCE_HIGH

    def test_gja_symbol_is_annotated_connexin(self):
        hit = cfa.SearchHit(
            database="ncbi_gene",
            level="A",
            query="q",
            gene_name="GJA1",
            protein_name="gap junction protein",
        )
        result = cfa.classify_hit("connexin", hit)
        assert result is not None
        assert result.final_category == cfa.CATEGORY_ANNOTATED
        assert result.matching_term == "GJA1"

    def test_gap_junction_is_related_not_annotated(self):
        hit = cfa.SearchHit(
            database="uniprot",
            level="B",
            query="q",
            protein_name="Gap junction protein homolog",
            gene_name="UNKNOWN1",
        )
        result = cfa.classify_hit("connexin", hit)
        assert result is not None
        assert result.final_category == cfa.CATEGORY_RELATED
        assert result.confidence == cfa.CONFIDENCE_MEDIUM

    def test_vague_channel_hit_ignored(self):
        hit = cfa.SearchHit(
            database="ncbi_protein",
            level="B",
            query="q",
            protein_name="Putative channel protein",
            gene_name="",
        )
        assert cfa.classify_hit("connexin", hit) is None

    def test_innexin_alias_shakb_annotated(self):
        hit = cfa.SearchHit(
            database="ncbi_gene",
            level="B",
            query="q",
            gene_name="shakB",
            protein_name="Protein shakB",
        )
        result = cfa.classify_hit("innexin", hit)
        assert result is not None
        assert result.final_category == cfa.CATEGORY_ANNOTATED

    def test_panx_symbol_annotated(self):
        hit = cfa.SearchHit(
            database="uniprot",
            level="A",
            query="q",
            gene_name="PANX1",
            protein_name="Pannexin-1",
        )
        result = cfa.classify_hit("pannexin", hit)
        assert result is not None
        assert result.final_category == cfa.CATEGORY_ANNOTATED

    def test_merge_prefers_annotated_over_related(self):
        annotated = cfa.Classification(
            final_category=cfa.CATEGORY_ANNOTATED,
            confidence=cfa.CONFIDENCE_HIGH,
            matching_term="connexin",
            database_source="uniprot",
            level="A",
        )
        related = cfa.Classification(
            final_category=cfa.CATEGORY_RELATED,
            confidence=cfa.CONFIDENCE_MEDIUM,
            matching_term="gap junction",
            database_source="ncbi_gene",
            level="B",
        )
        final = cfa.merge_family_results([annotated], [related])
        assert final.final_category == cfa.CATEGORY_ANNOTATED

    def test_merge_none_when_empty(self):
        final = cfa.merge_family_results([], [])
        assert final.final_category == cfa.CATEGORY_NONE


class TestHelpers:
    def test_gene_matches_prefix(self):
        assert cfa.gene_matches_prefix("GJA1", "GJA")
        assert not cfa.gene_matches_prefix("GJB2", "GJA")

    def test_dedupe_species(self):
        assert cfa.dedupe_preserve_order(["A", "a", "B"]) == ["A", "B"]

    def test_level_a_queries_include_symbol_patterns(self):
        queries = cfa.level_a_queries("connexin", "Homo sapiens")
        joined = " | ".join(q for _, q in queries)
        assert "gene:GJA*" in joined
        assert "protein_name:connexin" in joined

    def test_level_b_queries_include_gap_junction(self):
        queries = cfa.level_b_queries("innexin", "Drosophila melanogaster")
        joined = " | ".join(q for _, q in queries)
        assert "gap junction" in joined
        assert "shakB" in joined or "shakb" in joined.lower()

    def test_load_species_from_text(self, tmp_path):
        path = tmp_path / "species.txt"
        path.write_text("# comment\nHomo sapiens\nDanio rerio\n", encoding="utf-8")
        assert cfa.load_species_from_text(path) == ["Homo sapiens", "Danio rerio"]

    def test_write_results_and_summary(self, tmp_path):
        rows = [
            {
                "organism": "Homo sapiens",
                "family": "connexin",
                "exact_annotation_found": "yes",
                "related_annotation_found": "no",
                "final_category": cfa.CATEGORY_ANNOTATED,
                "confidence": "high",
                "top_matching_term": "connexin",
                "database_source": "uniprot",
                "search_level": "A",
                "top_accession": "P17302",
                "top_gene_symbol": "GJA1",
                "top_protein_name": "Gap junction alpha-1 protein",
                "level_a_hit_count": "3",
                "level_b_hit_count": "0",
                "example_query": "q1",
            },
            {
                "organism": "Homo sapiens",
                "family": "innexin",
                "exact_annotation_found": "no",
                "related_annotation_found": "no",
                "final_category": cfa.CATEGORY_NONE,
                "confidence": "",
                "top_matching_term": "",
                "database_source": "",
                "search_level": "",
                "top_accession": "",
                "top_gene_symbol": "",
                "top_protein_name": "",
                "level_a_hit_count": "0",
                "level_b_hit_count": "0",
                "example_query": "",
            },
        ]
        detail = tmp_path / "detail.csv"
        summary = tmp_path / "summary.csv"
        cfa.write_results(rows, detail)
        cfa.write_summary(rows, summary)
        assert "connexin" in detail.read_text(encoding="utf-8")
        assert "innexin_category" in summary.read_text(encoding="utf-8")


class TestSearchIntegration:
    def test_classify_species_family_mocked(self):
        annotated_hit = cfa.SearchHit(
            database="uniprot",
            level="A",
            query="q",
            accession="P17302",
            gene_name="GJA1",
            protein_name="Gap junction alpha-1 protein",
        )
        with patch.object(cfa, "run_level", side_effect=[
            ([cfa.classify_hit("connexin", annotated_hit)], 1),
            ([], 0),
        ]):
            row = cfa.classify_species_family("Homo sapiens", "connexin")
        assert row["final_category"] == cfa.CATEGORY_ANNOTATED
        assert row["exact_annotation_found"] == "yes"
