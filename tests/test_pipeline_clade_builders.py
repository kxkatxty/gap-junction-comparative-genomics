from __future__ import annotations

from pathlib import Path

import pipeline.innexin_clade_comparison_builder as inx
import pipeline.connexin_clade_comparison_builder as cnx


class TestInnexinCladeHelpers:
    def test_f_i(self):
        assert inx._f("") == 0.0
        assert inx._i("3.2") == 3

    def test_infer_clade(self):
        assert inx.infer_clade("Drosophila_melanogaster", "Drosophila melanogaster") == "Insecta"
        assert inx.infer_clade("Caenorhabditis_elegans", "Caenorhabditis elegans") == "Nematoda"
        assert inx.infer_clade("Adineta_vaga", "Adineta vaga") == "Rotifera"
        assert inx.infer_clade("Mytilus_edulis", "Mytilus edulis") == "Mollusca"
        assert inx.infer_clade("Mnemiopsis_leidyi", "Mnemiopsis leidyi") == "Ctenophora"
        assert inx.infer_clade("Unknown_thing", "Unknown thing") == "Other"

    def test_classify_reference_type(self):
        assert inx.classify_reference_type("shakB") == "shakB"
        assert inx.classify_reference_type("INX1_DROME") == "Inx1/ogre"
        assert inx.classify_reference_type("INX2_DROME") == "Inx2"
        assert inx.classify_reference_type("INX2_CAEEL") == "inx-2 (Cele)"
        assert inx.classify_reference_type("UNC-7_CAEEL") == "unc-7"
        assert inx.classify_reference_type("UNC9_CAEEL") == "unc-9"
        assert inx.classify_reference_type("something") == "other/unknown"

    def test_is_kept(self):
        assert inx._is_kept("rejected_false_positive") is False
        assert inx._is_kept("high_confidence_innexin_candidate") is True

    def test_key_findings_nonempty(self):
        rows = [
            {
                "species_dir": "Drosophila_melanogaster",
                "organism": "Drosophila melanogaster",
                "clade": "Insecta",
                "reference_type": "Inx2",
                "subfamily": "SF4_Inx2_expansion",
                "source": "reference_db",
                "protein_length": 367,
                "reference_identity": 1.0,
                "rank_category": "reference",
            },
            {
                "species_dir": "Caenorhabditis_elegans",
                "organism": "Caenorhabditis elegans",
                "clade": "Nematoda",
                "reference_type": "unc-9",
                "subfamily": "SF3_Inx3_Inx7_nematode",
                "source": "reference_db",
                "protein_length": 390,
                "reference_identity": 1.0,
                "rank_category": "reference",
            },
            {
                "species_dir": "Adineta_vaga",
                "organism": "Adineta vaga",
                "clade": "Rotifera",
                "reference_type": "unc-9",
                "subfamily": "SF3_Inx3_Inx7_nematode",
                "source": "curator_new",
                "protein_length": 340,
                "reference_identity": 0.4,
                "rank_category": "curator_present_innexin",
            },
        ]
        findings = inx.key_findings(rows)
        assert isinstance(findings, list)
        assert len(findings) >= 1


class TestConnexinCladeHelpers:
    def test_key_findings_nonempty(self):
        rows = [
            {
                "species_dir": "Homo_sapiens",
                "organism": "Homo sapiens",
                "clade": "Mammal",
                "reference_type": "GJA1",
                "subfamily": "alpha_GJA",
                "source": "reference_db",
                "protein_length": 382,
                "reference_identity": 1.0,
            },
            {
                "species_dir": "Gallus_gallus",
                "organism": "Gallus gallus",
                "clade": "Bird",
                "reference_type": "GJB2",
                "subfamily": "beta_GJB",
                "source": "reference_db",
                "protein_length": 226,
                "reference_identity": 1.0,
            },
        ]
        findings = cnx.key_findings(rows)
        assert isinstance(findings, list)
        assert len(findings) >= 1

    def test_plot_loci_by_clade_writes_png(self, tmp_path: Path):
        rows = [
            {"clade": "Mammal", "reference_type": "GJA1", "source": "reference_db"},
            {"clade": "Mammal", "reference_type": "GJB2", "source": "reference_db"},
            {"clade": "Bird", "reference_type": "GJA1", "source": "discovery"},
        ]
        out = tmp_path / "loci.png"
        cnx.plot_loci_by_clade(rows, out)
        assert out.exists()
        assert out.stat().st_size > 0
