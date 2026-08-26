from __future__ import annotations

from pathlib import Path

import numpy as np

import pipeline.innexin_similarity_builder as inx
import pipeline.connexin_similarity_builder as cnx


class TestInnexinSimilarityHelpers:
    def test_clean_aa_and_sf(self):
        assert inx._clean_aa("mk*xx") == "MK"
        assert inx._sf_from_gene("Inx2", "Drosophila melanogaster") == "SF4_Inx2_expansion"
        assert inx._sf_from_gene("inx-2", "Caenorhabditis elegans") == "SF3_Inx3_Inx7_nematode"
        assert inx._sf_from_gene("inx-7", "Caenorhabditis elegans") == "SF3_Inx3_Inx7_nematode"

    def test_read_fasta(self, tmp_path: Path):
        path = tmp_path / "a.fa"
        path.write_text(">one\nAAA\n>two\nCCC\n", encoding="utf-8")
        assert inx._read_fasta(path) == [("one", "AAA"), ("two", "CCC")]
        assert inx._read_fasta(tmp_path / "missing.fa") == []

    def test_cosine_and_aa_matrix(self):
        rows = [
            {"seq_id": "a", "seq": "ACDEFGHIKLMNPQRSTVWY" * 2},
            {"seq_id": "b", "seq": "ACDEFGHIKLMNPQRSTVWY" * 2},
            {"seq_id": "c", "seq": "WWWWWWWWWWWWWWWWWWWW" * 2},
        ]
        X, ids = inx.aa_composition_matrix(rows)
        assert ids == ["a", "b", "c"]
        sim = inx.cosine_sim(X)
        assert sim.shape == (3, 3)
        assert sim[0, 1] > sim[0, 2]

    def test_kmer_embedding(self):
        rows = [
            {"seq_id": "a", "seq": "ACDEFGHIKLMNPQRSTVWY" * 5},
            {"seq_id": "b", "seq": "ACDEFGHIKLMNPQRSTVWY" * 5},
        ]
        emb, ids = inx.kmer_embedding_matrix(rows, k=3, n_comp=4)
        assert len(ids) == 2
        assert emb.shape[0] == 2

    def test_summarize_clusters_and_pairwise(self, tmp_path: Path, monkeypatch):
        monkeypatch.setattr(inx, "OUT", tmp_path)
        clu = tmp_path / "clu.tsv"
        clu.write_text("rep1\tmem1\nrep1\tmem2\nrep2\tmem3\n", encoding="utf-8")
        meta = {
            "mem1": {"species_dir": "SpA", "subfamily": "SF1"},
            "mem2": {"species_dir": "SpB", "subfamily": "SF1"},
            "mem3": {"species_dir": "SpC", "subfamily": "SF2"},
            "rep1": {"species_dir": "SpA", "subfamily": "SF1"},
            "rep2": {"species_dir": "SpC", "subfamily": "SF2"},
        }
        summary = inx.summarize_clusters({0.5: clu}, meta)
        assert summary[0]["n_clusters"] == 2
        assert summary[0]["n_singletons"] == 1

        m8 = tmp_path / "hits.m8"
        m8.write_text(
            "mem1\tmem2\t80.0\t100\nmem1\tmem3\t30.0\t100\n",
            encoding="utf-8",
        )
        pw = inx.pairwise_from_m8(m8, meta)
        assert len(pw["within"]["SF1"]) == 1
        assert pw["within"]["SF1"][0] == 80.0


class TestConnexinSimilarityHelpers:
    def test_aa_and_kmer_matrix(self):
        rows = [
            {"seq_id": "a", "seq": "ACDEFGHIKLMNPQRSTVWY" * 3, "subfamily": "alpha_GJA", "species_dir": "Hs"},
            {"seq_id": "b", "seq": "ACDEFGHIKLMNPQRSTVWY" * 3, "subfamily": "alpha_GJA", "species_dir": "Mm"},
        ]
        X, ids = cnx.aa_matrix(rows)
        assert len(ids) == 2
        emb, ids2 = cnx.kmer_matrix(rows, k=3, n_comp=3)
        assert len(ids2) == 2
        assert emb.shape[0] == 2

    def test_group_table(self):
        rows = [
            {"seq_id": "a", "seq": "A" * 20, "subfamily": "alpha_GJA", "species_dir": "Hs"},
            {"seq_id": "b", "seq": "A" * 20, "subfamily": "alpha_GJA", "species_dir": "Mm"},
            {"seq_id": "c", "seq": "W" * 20, "subfamily": "beta_GJB", "species_dir": "Gg"},
        ]
        X, ids = cnx.aa_matrix(rows)
        # simple identity-like sim
        sim = X @ X.T
        meta = {r["seq_id"]: r for r in rows}
        table = cnx.group_table(sim, ids, meta, "aa")
        assert any(r["same_group"] == "true" for r in table)

    def test_findings(self):
        sweep = [{"min_seq_id": 0.5, "n_clusters": 2, "mean_subfamily_purity": 0.9, "n_multi_clusters": 1}]
        pw = {"within": {"alpha_GJA": [70.0]}, "between": {("ANY_BETWEEN",): [30.0]}, "pairs": []}
        finds = cnx.findings(sweep, pw, [])
        assert isinstance(finds, list)
        assert len(finds) >= 1
