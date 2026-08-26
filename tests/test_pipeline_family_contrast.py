from __future__ import annotations

from pathlib import Path

import numpy as np
import pandas as pd

import pipeline.family_comparison_builder as fam
import pipeline.inx_vs_cnx_builder as ivc


class TestFamilyComparisonHelpers:
    def test_boolish_and_short(self):
        assert fam._boolish("true") is True
        assert fam._boolish("YES") is True
        assert fam._boolish("0") is False
        assert fam._short("SF3_Inx3_Inx7_nematode", "innexin")
        assert fam._short("alpha_GJA", "connexin")

    def test_weighted_mean_and_metrics(self):
        df = pd.DataFrame(
            [
                {"method": "mmseqs", "group_a": "A", "group_b": "A", "same_group": True, "mean_similarity": 50.0, "median_similarity": 48.0, "n_pairs": 10},
                {"method": "mmseqs", "group_a": "A", "group_b": "B", "same_group": False, "mean_similarity": 20.0, "median_similarity": 18.0, "n_pairs": 5},
                {"method": "mmseqs", "group_a": "B", "group_b": "B", "same_group": True, "mean_similarity": 60.0, "median_similarity": 55.0, "n_pairs": 10},
            ]
        )
        metrics = fam.family_metrics(df, "mmseqs")
        assert metrics["n_within_pairs"] == 20
        assert metrics["n_between_pairs"] == 5
        assert metrics["within_mean"] > metrics["between_mean"]

    def test_pair_mean_and_heatmap(self):
        df = pd.DataFrame(
            [
                {"method": "mmseqs", "group_a": "A", "group_b": "B", "mean_similarity": 33.0, "n_pairs": 4, "same_group": False, "median_similarity": 30.0},
                {"method": "mmseqs", "group_a": "A", "group_b": "A", "mean_similarity": 70.0, "n_pairs": 2, "same_group": True, "median_similarity": 70.0},
            ]
        )
        assert fam.pair_mean(df, "mmseqs", "A", "B") == 33.0
        mat = fam.heatmap_matrix(df, "mmseqs", ["A", "B"])
        assert mat.shape == (2, 2)
        assert mat[0, 1] == 33.0

    def test_load_group_matrix(self, tmp_path: Path):
        path = tmp_path / "m.csv"
        path.write_text(
            "method,group_a,group_b,same_group,mean_similarity,median_similarity,n_pairs\n"
            "mmseqs,A,A,true,50,49,3\n",
            encoding="utf-8",
        )
        df = fam.load_group_matrix(path)
        assert bool(df.iloc[0]["same_group"]) is True
        assert df.iloc[0]["n_pairs"] == 3


class TestInxVsCnxHelpers:
    def test_seq_row_and_read_fasta(self, tmp_path: Path):
        row = ivc._seq_row("s1", "ACDEFGHIKLMNPQRSTVWYC", "innexin")
        assert row["length"] == 21
        assert row["family"] == "innexin"
        assert row["cys_per_100"] > 0
        assert "aa_A" in row

        fasta = tmp_path / "x.fa"
        fasta.write_text(">a\nACDE\n>b note\nFGHI\n", encoding="utf-8")
        rows = ivc.read_fasta(fasta, "connexin")
        assert [r["seq_id"] for r in rows] == ["a", "b"]
        assert rows[0]["family"] == "connexin"

    def test_kingdom(self):
        assert ivc.kingdom("Mammal") == "Vertebrates"
        assert ivc.kingdom("Ray-finned fish") == "Vertebrates"
        assert ivc.kingdom("Insecta") == "Invertebrates"

    def test_aa_matrix_pca_cosine(self):
        rows = [
            ivc._seq_row("a", "ACDEFGHIKLMNPQRSTVWY" * 3, "innexin"),
            ivc._seq_row("b", "ACDEFGHIKLMNPQRSTVWY" * 3, "innexin"),
            ivc._seq_row("c", "ACDEFGHIKLMNPQRSTVWY" * 3, "innexin"),
            ivc._seq_row("d", "WWWWWWWWWWWWWWWWWWWW" * 3, "connexin"),
            ivc._seq_row("e", "WWWWWWWWWWWWWWWWWWWW" * 3, "connexin"),
        ]
        X = ivc.aa_matrix(rows)
        assert X.shape[0] == 5
        xy = ivc.pca2(X)
        assert xy.shape == (5, 2)
        same = ivc.mean_cosine(X, np.array([0, 1, 2]), np.array([0, 1, 2]))
        cross = ivc.mean_cosine(X, np.array([0, 1, 2]), np.array([3, 4]))
        assert not np.isnan(same)
        assert not np.isnan(cross)
        assert same > cross
        assert np.isnan(ivc.mean_cosine(X, np.array([]), np.array([0])))

    def test_kmer_svd_shape(self):
        seqs = ["ACDEFGHIKLMNPQRSTVWY" * 4, "ACDEFGHIKLMNPQRSTVWY" * 4, "WWWWWWWWWWWWWWWWWWWW" * 4]
        emb = ivc.kmer_svd(seqs, k=3, n_comp=4)
        assert emb.shape[0] == 3
        assert emb.shape[1] <= 4
