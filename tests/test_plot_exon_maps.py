from __future__ import annotations

from pathlib import Path

import pandas as pd
import plot_exon_maps as pem


class TestParsing:
    def test_parse_attributes(self):
        attrs = pem.parse_attributes("Parent=tx1;phase=0")
        assert attrs["Parent"] == "tx1"

    def test_read_fasta_sequence(self, sample_fasta):
        assert pem.read_fasta_sequence(sample_fasta) == "MKTAYIAKQRQISFVKSHFSRQ"

    def test_parse_feature_line_exon(self):
        line = "chr1\t.\texon\t200\t400\t.\t+\t.\tParent=tx1\n"
        parsed = pem._parse_feature_line(line, "tx1")
        assert parsed == ("exon", (200, 400))

    def test_parse_feature_line_cds(self):
        line = "chr1\t.\tcds\t200\t400\t.\t+\t1\tParent=tx1\n"
        parsed = pem._parse_feature_line(line, "tx1")
        assert parsed == ("cds", (200, 400, 1))

    def test_parse_feature_line_wrong_parent(self):
        line = "chr1\t.\tcds\t200\t400\t.\t+\t0\tParent=tx2\n"
        assert pem._parse_feature_line(line, "tx1") is None


class TestCdsMapping:
    def test_map_cds_to_aa_plus_strand(self):
        cds = [(200, 400, 0), (800, 900, 0)]
        pieces, length = pem.map_cds_to_aa(cds, "+")
        assert length == 67 + 33
        assert pieces[0] == (200, 400, 1, 67)
        assert pieces[1] == (800, 898, 68, 100)

    def test_map_cds_to_aa_with_phase(self):
        cds = [(100, 199, 2)]
        _, length = pem.map_cds_to_aa(cds, "+")
        segment_length = 199 - 100 + 1
        assert length == (segment_length - 2) // 3

    def test_map_cds_to_aa_minus_strand(self):
        cds = [(800, 900, 0), (200, 400, 0)]
        pieces, length = pem.map_cds_to_aa(cds, "-")
        assert length > 0
        assert pieces[0][2] == 1

    def test_genomic_overlap_to_aa(self):
        aa_start, aa_end = pem.genomic_overlap_to_aa(210, 230, 200, 1)
        assert aa_start == 1 + (210 - 200) // 3
        assert aa_end == 1 + (230 - 200) // 3

    def test_coding_exon_aa_ranges(self, mini_gff3):
        exons, cds = pem.parse_transcript_cds_and_exons(mini_gff3, "tx1")
        ranges = pem.coding_exon_aa_ranges(exons, cds, "+", protein_length=100)
        assert len(ranges) == 2
        assert ranges[0]["coding_exon_number"] == 1
        assert ranges[0]["aa_start"] == 1
        assert ranges[1]["aa_end"] == 100


class TestTranscriptSelection:
    def test_pick_representative_transcripts(self):
        df = pd.DataFrame(
            [
                {
                    "family": "connexin",
                    "organism": "A",
                    "fasta_path": "a.fasta",
                    "uniprot_accession": "P1",
                    "transcript_id": "tx_long",
                    "reference_gene_symbol": "GJA1",
                    "gene_symbol": "GJA1",
                    "exon_number": 1,
                    "exon_start": 1,
                    "exon_end": 100,
                },
                {
                    "family": "connexin",
                    "organism": "A",
                    "fasta_path": "a.fasta",
                    "uniprot_accession": "P1",
                    "transcript_id": "tx_long",
                    "reference_gene_symbol": "GJA1",
                    "gene_symbol": "GJA1",
                    "exon_number": 2,
                    "exon_start": 200,
                    "exon_end": 400,
                },
                {
                    "family": "connexin",
                    "organism": "A",
                    "fasta_path": "a.fasta",
                    "uniprot_accession": "P1",
                    "transcript_id": "tx_short",
                    "reference_gene_symbol": "GJA1",
                    "gene_symbol": "GJA1",
                    "exon_number": 1,
                    "exon_start": 1,
                    "exon_end": 50,
                },
            ]
        )
        picked = pem.pick_representative_transcripts(df)
        assert set(picked["transcript_id"]) == {"tx_long"}

    def test_chunk_keys(self):
        keys = [(f"org{i}", "path", f"tx{i}", "GJA1", "GJA1") for i in range(5)]
        chunks = pem.chunk_keys(keys, 2)
        assert len(chunks) == 3
        assert len(chunks[0]) == 2


class TestFeatureCombining:
    def test_combine_prefers_seth_disorder(self):
        classic = pd.DataFrame(
            [
                {
                    "uniprot_accession": "P17302",
                    "track": "disorder",
                    "feature_id": "disorder",
                    "aa_start": 1,
                    "aa_end": 10,
                },
                {
                    "uniprot_accession": "P17302",
                    "track": "pfam",
                    "feature_id": "PF00822",
                    "aa_start": 1,
                    "aa_end": 100,
                },
            ]
        )
        rostlab = pd.DataFrame(
            [
                {
                    "uniprot_accession": "P17302",
                    "track": "disorder",
                    "feature_id": "disorder",
                    "aa_start": 50,
                    "aa_end": 80,
                    "source": "biocentral_seth",
                },
                {
                    "uniprot_accession": "P17302",
                    "track": "transmembrane",
                    "feature_id": "tm_helix",
                    "aa_start": 20,
                    "aa_end": 40,
                },
            ]
        )
        combined, has_rostlab = pem.combine_protein_features(classic, rostlab, "P17302")
        assert has_rostlab
        tracks = set(combined["track"])
        assert "pfam" in tracks
        assert "transmembrane" in tracks
        disorder = combined[combined["track"] == "disorder"]
        assert len(disorder) == 1
        assert disorder.iloc[0]["aa_start"] == 50


class TestLabels:
    def test_label_for_row(self):
        row = pd.Series(
            {
                "reference_gene_symbol": "GJA1",
                "gene_symbol": "gja1",
                "organism": "Danio rerio",
            }
        )
        assert "Danio rerio" in pem.label_for_row(row)

    def test_short_organism(self):
        assert pem.short_organism("Homo sapiens") == "H. sapiens"
        assert pem.short_organism("Drosophila") == "Drosophila"

    def test_brighten_color(self):
        assert pem._brighten_color("#E67E22") == "#FFB020"
        assert pem._brighten_color("#ABCDEF") == "#ABCDEF"


class TestPlotOutput:
    def test_plot_exon_count_boxplot_writes_file(self, tmp_path: Path):
        df = pd.DataFrame(
            {
                "family": ["connexin", "connexin", "innexin"],
                "fasta_path": ["a.fasta", "b.fasta", "c.fasta"],
                "exon_count": [3, 5, 2],
            }
        )
        out = tmp_path / "boxplot.png"
        pem.plot_exon_count_boxplot(df, out)
        assert out.exists()
        assert out.stat().st_size > 0

    def test_plot_compressed_exon_map_writes_file(self, tmp_path: Path):
        df = pd.DataFrame(
            [
                {
                    "family": "connexin",
                    "organism": "Homo sapiens",
                    "fasta_path": "a.fasta",
                    "uniprot_accession": "P17302",
                    "transcript_id": "tx1",
                    "reference_gene_symbol": "GJA1",
                    "gene_symbol": "GJA1",
                    "exon_number": 1,
                    "exon_start": 1,
                    "exon_end": 100,
                    "exon_length": 100,
                    "strand": "+",
                },
                {
                    "family": "connexin",
                    "organism": "Homo sapiens",
                    "fasta_path": "a.fasta",
                    "uniprot_accession": "P17302",
                    "transcript_id": "tx1",
                    "reference_gene_symbol": "GJA1",
                    "gene_symbol": "GJA1",
                    "exon_number": 2,
                    "exon_start": 200,
                    "exon_end": 300,
                    "exon_length": 101,
                    "strand": "+",
                },
            ]
        )
        out = tmp_path / "compressed.png"
        pem.plot_compressed_exon_map(df, "connexin", out)
        assert out.exists()

    def test_plot_genomic_exon_map_writes_file(self, tmp_path: Path):
        df = pd.DataFrame(
            [
                {
                    "family": "connexin",
                    "organism": "Homo sapiens",
                    "fasta_path": "a.fasta",
                    "uniprot_accession": "P17302",
                    "transcript_id": "tx1",
                    "reference_gene_symbol": "GJA1",
                    "gene_symbol": "GJA1",
                    "exon_number": 1,
                    "exon_start": 100,
                    "exon_end": 200,
                    "exon_length": 101,
                    "strand": "+",
                }
            ]
        )
        out = tmp_path / "genomic.png"
        pem.plot_genomic_exon_map(df, "connexin", out)
        assert out.exists()

    def test_normalize_species_name(self):
        assert pem.normalize_species_name("Homo sapiens") == "Homo_sapiens"

    def test_transcript_keys_sorted(self):
        df = pd.DataFrame(
            {
                "organism": ["B", "A"],
                "fasta_path": ["b.fasta", "a.fasta"],
                "transcript_id": ["tx2", "tx1"],
                "reference_gene_symbol": ["GJA1", "GJA1"],
                "gene_symbol": ["GJA1", "GJA1"],
            }
        )
        keys = pem.transcript_keys(df)
        assert keys[0][0] == "A"

