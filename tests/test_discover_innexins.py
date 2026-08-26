from __future__ import annotations

from pathlib import Path

import discover_innexins as di


class TestFastaAndSpecies:
    def test_read_species_txt(self, tmp_path: Path):
        path = tmp_path / "species.txt"
        path.write_text("# comment\nPriapulus caudatus\n\nAdineta vaga\n", encoding="utf-8")
        assert di.read_species_txt(path) == ["Priapulus caudatus", "Adineta vaga"]

    def test_read_and_write_fasta(self, tmp_path: Path):
        records = {"seqA": "ACGTACGT", "seqB": "TTTT"}
        path = tmp_path / "test.fasta"
        di.write_fasta_records(records, path)
        assert di.read_fasta_records(path) == records


class TestClustering:
    def test_cluster_hsps_same_scaffold(self):
        hsps = [
            di.BlastHsp("q1", "scaf1", 1, 100, 1000, 1800, 1e-20, 200, 45, 1),
            di.BlastHsp("q2", "scaf1", 1, 100, 5000, 5600, 1e-15, 150, 40, 1),
            di.BlastHsp("q3", "scaf1", 1, 100, 20000, 20800, 1e-10, 120, 35, 1),
        ]
        loci = di.cluster_hsps_into_loci(hsps, max_gap_bp=10000)
        assert len(loci) == 2

    def test_hsp_quality_filters_reject_short_or_weak(self):
        weak = di.BlastHsp("q1", "scaf1", 1, 15, 1000, 1080, 1e-12, 80, 40, 1)
        short_bp = di.BlastHsp("q1", "scaf1", 1, 100, 1000, 1050, 1e-12, 80, 40, 1)
        low_bits = di.BlastHsp("q1", "scaf1", 1, 100, 1000, 1200, 1e-12, 20, 40, 1)
        low_pident = di.BlastHsp("q1", "scaf1", 1, 100, 1000, 1200, 1e-12, 80, 15, 1)
        short_low_pident = di.BlastHsp("q1", "scaf1", 1, 30, 1000, 1200, 1e-12, 80, 20, 1)
        good = di.BlastHsp("q1", "scaf1", 1, 100, 1000, 1200, 1e-12, 80, 40, 1)
        assert not di.passes_hsp_quality_filters(weak)
        assert not di.passes_hsp_quality_filters(short_bp)
        assert not di.passes_hsp_quality_filters(low_bits)
        assert not di.passes_hsp_quality_filters(low_pident)
        assert not di.passes_hsp_quality_filters(short_low_pident)
        assert di.passes_hsp_quality_filters(good)

    def test_build_exon_model_merges_close_hsps(self):
        hsps = [
            di.BlastHsp("q1", "scaf1", 1, 80, 100, 400, 1e-10, 100, 40, 1),
            di.BlastHsp("q1", "scaf1", 90, 150, 900, 1200, 1e-10, 90, 38, 1),
        ]
        exons = di.build_exon_model_from_hsps(hsps)
        assert len(exons) == 2
        assert exons[0].start == 100


class TestTranslationAndTM:
    def test_translate_dna_start_methionine(self):
        # ATG AAA TAA
        protein = di.translate_dna("ATGAAATAA", frame=1)
        assert protein == "MK"

    def test_predict_tm_helices_innexin_like(self):
        protein = "M" + "LALALALALALALALALAL" * 4 + "G" * 80
        assert di.predict_tm_helices(protein) >= 1


class TestValidation:
    def test_false_positive_keyword_detected(self):
        assert di.contains_false_positive_keyword("receptor tyrosine kinase") == "kinase"

    def test_score_high_confidence_candidate(self):
        cand = di.CandidateLocus(
            organism="Test",
            candidate_id="candidate_inx_001",
            seqid="scaf1",
            strand="+",
            locus_start=1,
            locus_end=2000,
            window_start=1,
            window_end=4000,
            protein_sequence="M" + "A" * 380,
            protein_length=381,
            exon_count=2,
            tm_helix_count=4,
            reference_identity=35.0,
            reference_coverage=70.0,
            reference_hit_count=3,
            completeness="complete",
        )
        score, category, _ = di.score_candidate(cand)
        assert score >= 60
        assert category == di.RANK_HIGH

    def test_score_rejects_short_protein(self):
        cand = di.CandidateLocus(
            organism="Test",
            candidate_id="candidate_inx_002",
            seqid="scaf1",
            strand="+",
            locus_start=1,
            locus_end=200,
            window_start=1,
            window_end=400,
            protein_sequence="M" + "A" * 50,
            protein_length=51,
            exon_count=1,
            tm_helix_count=0,
        )
        _, category, reason = di.score_candidate(cand)
        assert category == di.RANK_REJECTED
        assert reason


class TestMiniprotParsing:
    def test_parse_miniprot_gff(self, tmp_path: Path):
        path = tmp_path / "miniprot.gff"
        path.write_text(
            "scaf1\tminiprot\tmRNA\t1000\t1300\t584\t+\t.\t"
            "ID=MP000001;Rank=1;Identity=0.3511;Target=ref1 10 80\n"
            "scaf2\tminiprot\tmRNA\t4700\t5000\t310\t-\t.\t"
            "ID=MP000002;Rank=1;Identity=0.2697;Target=ref2 5 60\n",
            encoding="utf-8",
        )
        hsps = di.parse_miniprot_gff(path)
        assert len(hsps) == 2
        assert hsps[0].target_id == "scaf1"
        assert hsps[0].target_start == 1000
        assert hsps[0].identity == 35.11
        assert hsps[0].frame == 1
        assert hsps[1].frame == -1
        assert hsps[1].target_start == 4700
        assert hsps[1].target_end == 5000


class TestOutputs:
    def test_write_gff3_and_fasta(self, tmp_path: Path):
        cand = di.CandidateLocus(
            organism="Priapulus caudatus",
            candidate_id="candidate_inx_001",
            seqid="scaffold1",
            strand="+",
            locus_start=1000,
            locus_end=5000,
            window_start=800,
            window_end=5200,
            exons=[di.ExonModel(1, 1000, 1500), di.ExonModel(2, 3000, 5000)],
            protein_sequence="MKTAYIAK",
            protein_length=8,
            exon_count=2,
            tm_helix_count=4,
            rank_category=di.RANK_HIGH,
            rank_score=82.0,
        )
        gff = tmp_path / "out.gff3"
        fasta = tmp_path / "out.fasta"
        di.write_candidate_gff3([cand], gff)
        di.write_candidate_fasta([cand], fasta)
        gff_text = gff.read_text(encoding="utf-8")
        assert "candidate_inx_001" in gff_text
        assert "innexin_discovery" in gff_text
        assert "MKTAYIAK" in fasta.read_text(encoding="utf-8")
