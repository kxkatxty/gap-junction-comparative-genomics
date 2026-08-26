from __future__ import annotations

from pathlib import Path

import pipeline.connexin_common as cc


class TestConnexinCommonHelpers:
    def test_f_and_i_defaults(self):
        assert cc._f(None) == 0.0
        assert cc._f("") == 0.0
        assert cc._f("bad", default=3.5) == 3.5
        assert cc._f("2.5") == 2.5
        assert cc._i(None) == 0
        assert cc._i("4.9") == 4
        assert cc._i("x", default=7) == 7

    def test_clean_aa(self):
        assert cc.clean_aa("mk*TA") == "MK"
        assert cc.clean_aa("acgt!") == "ACGT"
        assert "X" in cc.clean_aa("ABJ")  # B/J not standard AA

    def test_classify_connexin_type(self):
        assert cc.classify_connexin_type("GJA1") == "GJA1"
        assert cc.classify_connexin_type("CX43") == "GJA1"
        assert cc.classify_connexin_type("GJB2_HUMAN") == "GJB2"
        assert cc.classify_connexin_type("GJC1") == "GJC"
        assert cc.classify_connexin_type("GJD2") == "GJD"
        assert cc.classify_connexin_type("GJE1") == "GJE"
        assert cc.classify_connexin_type("random") == "other/unknown"

    def test_infer_clade_from_name_hints(self, monkeypatch):
        monkeypatch.setattr(cc, "catalog", lambda: {})
        assert cc.infer_clade("Homo_sapiens", "Homo sapiens") == "Mammal"
        assert cc.infer_clade("Gallus_gallus", "Gallus gallus") == "Bird"
        assert cc.infer_clade("Danio_rerio", "Danio rerio") == "Ray-finned fish"
        assert cc.infer_clade("Xenopus_laevis", "Xenopus laevis") == "Amphibian"
        assert cc.infer_clade("Leucoraja_erinaceus", "Leucoraja erinaceus") == "Cartilaginous fish"
        assert cc.infer_clade("Mystery_sp", "Mystery sp") == "Other vertebrate"

    def test_infer_clade_uses_catalog(self, monkeypatch):
        monkeypatch.setattr(cc, "catalog", lambda: {"Homo sapiens": "Mammal", "Homo_sapiens": "Mammal"})
        assert cc.infer_clade("Homo_sapiens", "Homo sapiens") == "Mammal"

    def test_is_kept_and_high_conf(self):
        assert cc.is_kept_rank("high_confidence_connexin_candidate") is True
        assert cc.is_kept_rank("weak_manual_review") is True
        assert cc.is_kept_rank("rejected_false_positive") is False
        assert cc.is_high_conf("high_confidence_connexin_candidate") is True
        assert cc.is_high_conf("reference") is True
        assert cc.is_high_conf("weak_manual_review") is False

    def test_read_write_fasta(self, tmp_path: Path):
        fasta = tmp_path / "in.fa"
        fasta.write_text(">seq1\nACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWYACDEFGHIKLMNPQRSTVWY\n>seq2\nSHORT\n", encoding="utf-8")
        records = cc.read_fasta(fasta)
        assert len(records) == 2
        assert records[0][0] == "seq1"

        out = tmp_path / "out.fa"
        rows = [
            {"seq_id": "ok", "seq": "ACDEFGHIKLMNPQRSTVWY" * 5},
            {"seq_id": "too_short", "seq": "ACDE"},
        ]
        cc.write_fasta(out, rows)
        text = out.read_text(encoding="utf-8")
        assert ">ok" in text
        assert "too_short" not in text

    def test_apply_best_hits(self, tmp_path: Path):
        m8 = tmp_path / "hits.m8"
        m8.write_text(
            "disc__Homo_sapiens__c1\tref__Homo_sapiens__GJA1__P17302\t85.0\t100\t0\t0\t1\t100\t1\t100\t1e-40\t200\n",
            encoding="utf-8",
        )
        rows = [{"seq_id": "disc__Homo_sapiens__c1", "reference_type": "other/unknown", "subfamily": "unassigned", "gene_label": ""}]
        cc.apply_best_hits(rows, m8)
        assert rows[0]["reference_type"] == "GJA1"
        assert rows[0]["gene_label"] == "GJA1"
        assert rows[0]["reference_identity"] == 0.85

    def test_summarize(self):
        rows = [
            {"species_dir": "A", "clade": "Mammal", "reference_type": "GJA1", "source": "reference_db", "subfamily": "alpha_GJA"},
            {"species_dir": "A", "clade": "Mammal", "reference_type": "GJB2", "source": "discovery", "subfamily": "beta_GJB"},
            {"species_dir": "B", "clade": "Bird", "reference_type": "GJA1", "source": "discovery", "subfamily": "alpha_GJA"},
        ]
        stats = cc.summarize(rows)
        assert stats["n_loci"] == 3
        assert stats["n_species"] == 2
        assert stats["by_clade"]["Mammal"] == 2
        assert stats["by_type"]["GJA1"] == 2
