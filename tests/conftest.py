from __future__ import annotations

import matplotlib

matplotlib.use("Agg")

from pathlib import Path

import pytest


FIXTURES = Path(__file__).parent / "fixtures"


@pytest.fixture
def mini_gff3(tmp_path: Path) -> Path:
    content = """##gff-version 3
chr1\tRefSeq\tgene\t100\t5000\t.\t+\t.\tID=gene1;Name=GJA1;gene_biotype=protein_coding;product=gap junction alpha-1 protein
chr1\tRefSeq\tmRNA\t200\t4900\t.\t+\t.\tID=tx1;Parent=gene1;product=gap junction alpha-1 protein
chr1\tRefSeq\texon\t200\t400\t.\t+\t.\tParent=tx1
chr1\tRefSeq\texon\t800\t1200\t.\t+\t.\tParent=tx1
chr1\tRefSeq\texon\t3000\t4900\t.\t+\t.\tParent=tx1
chr1\tRefSeq\tcds\t200\t400\t.\t+\t0\tParent=tx1
chr1\tRefSeq\tcds\t800\t900\t.\t+\t0\tParent=tx1
chr1\tRefSeq\tgene\t6000\t8000\t.\t+\t.\tID=gene2;Name=cnst1;gene_biotype=protein_coding
chr1\tRefSeq\tmRNA\t6100\t7900\t.\t+\t.\tID=tx2;Parent=gene2
chr1\tRefSeq\texon\t6100\t7900\t.\t+\t.\tParent=tx2
chr1\tRefSeq\tgene\t9000\t10000\t.\t+\t.\tID=gene3;Name=GJB2;gene_biotype=protein_coding;product=gap junction beta-2 protein
chr1\tRefSeq\tmRNA\t9100\t9900\t.\t+\t.\tID=tx3;Parent=gene3
chr1\tRefSeq\texon\t9100\t9900\t.\t+\t.\tParent=tx3
"""
    path = tmp_path / "mini.gff3"
    path.write_text(content, encoding="utf-8")
    return path


@pytest.fixture
def sample_fasta(tmp_path: Path) -> Path:
    content = (
        ">sp|P17302|CXA1_HUMAN Gap junction alpha-1 protein OS=Homo sapiens GN=GJA1\n"
        "MKTAYIAKQRQISFVKSHFSRQ\n"
    )
    path = tmp_path / "GJA1__P17302.fasta"
    path.write_text(content, encoding="utf-8")
    return path


@pytest.fixture
def whitelist_csv(tmp_path: Path) -> Path:
    content = """plot_group,family,reference_gene_symbol,uniprot_accession,organism,enabled,title_label,notes
GJA1_human,connexin,GJA1,P17302,Homo sapiens,true,GJA1,Human only
GJA1,connexin,GJA1,,,true,GJA1,All species
"""
    path = tmp_path / "gene_plot_whitelist.csv"
    path.write_text(content, encoding="utf-8")
    return path


@pytest.fixture
def mapping_csv(tmp_path: Path) -> Path:
    content = """family,organism,fasta_path,uniprot_accession,reference_gene_symbol
connexin,Homo sapiens,project/data/references/connexins/Homo_sapiens/GJA1__P17302.fasta,P17302,GJA1
connexin,Danio rerio,project/data/references/connexins/Danio_rerio/gja1__P18246.fasta,P18246,gja1
"""
    path = tmp_path / "protein_mapping.csv"
    path.write_text(content, encoding="utf-8")
    return path


@pytest.fixture
def pipeline_sandbox(tmp_path: Path, monkeypatch):
    """Minimal project tree for pipeline module tests."""
    root = tmp_path
    metadata = root / "project" / "metadata"
    results = root / "project" / "results"
    refs = root / "project" / "data" / "references"
    synvoy = root / "SynVoy"
    discovery = results / "innexin_discovery"
    for path in (
        metadata,
        results,
        refs / "innexins" / "Drosophila_melanogaster",
        discovery / "Test_species",
        synvoy / "results" / "innexin_synvoy" / "dmel_ogre",
        results / "plots",
    ):
        path.mkdir(parents=True, exist_ok=True)

    (metadata / "species_config.csv").write_text(
        "family,organism,reviewed_only\n"
        "innexin,Drosophila melanogaster,true\n"
        "connexin,Homo sapiens,true\n",
        encoding="utf-8",
    )
    (metadata / "not_annotated_innexin.txt").write_text("Test species\n", encoding="utf-8")
    (metadata / "family_annotation_status_test.csv").write_text(
        "organism,family,final_category,top_matching_term,final_evidence_level\n"
        "Test species,innexin,no_related_entry_found,,\n"
        "Drosophila melanogaster,innexin,annotated_family_member_found,Inx2,high\n",
        encoding="utf-8",
    )
    (metadata / "gap_junction_quality_table.csv").write_text(
        "organism,assembly_level,annotation_available\n"
        "Drosophila melanogaster,Chromosome,yes\n",
        encoding="utf-8",
    )
    (refs / "innexins" / "Drosophila_melanogaster" / "Inx2__Q9V427.fasta").write_text(
        ">sp|Q9V427|INX2_DROME Innexin inx2\n"
        "MKTAYIAKQRQISFVKSHFSRQ\n",
        encoding="utf-8",
    )
    (discovery / "Test_species" / "candidate_loci.csv").write_text(
        "organism,candidate_id,seqid,strand,locus_start,locus_end,exon_count,protein_length,"
        "tm_helix_count,reference_identity,reference_coverage,best_reference_hit,"
        "rank_category,rank_score,validation_flags\n"
        "Test species,candidate_inx_001,scaf1,+,100,500,2,350,4,45.0,60.0,Inx2,"
        "high_confidence_innexin_candidate,65.0,tm_helices_ok\n"
        "Test species,candidate_inx_002,scaf2,-,600,900,1,120,1,10.0,20.0,,"
        "rejected_false_positive,0.0,too_short\n",
        encoding="utf-8",
    )
    (discovery / "Test_species" / "candidate_innexins.fasta").write_text(
        ">candidate_inx_001 Test species\n"
        "MKTAYIAKQRQISFVKSHFSRQ\n",
        encoding="utf-8",
    )
    (discovery / "Test_species" / "discovery_summary.json").write_text(
        '{"organism": "Test species", "accepted": 1}',
        encoding="utf-8",
    )
    (results / "plots" / "innexin_compressed_exon_map.png").write_bytes(b"\x89PNG\r\n\x1a\n")
    (synvoy / "results" / "innexin_synvoy" / "dmel_ogre" / "synteny_block_1_synteny_plot.html").write_text(
        "<html></html>",
        encoding="utf-8",
    )

    import pipeline.candidate_validator as validator
    import pipeline.common as common
    import pipeline.family_result_merger as merger
    import pipeline.figure_export_pipeline as figures
    import pipeline.phylogeny_input_builder as phylo_input
    import pipeline.species_panel_builder as panel
    import pipeline.summary_table_builder as summary
    import pipeline.tree_annotation_helper as tree_helper

    discovery_dirs = (discovery,)
    monkeypatch.setattr(common, "PROJECT_ROOT", root)
    monkeypatch.setattr(common, "METADATA_DIR", metadata)
    monkeypatch.setattr(common, "RESULTS_DIR", results)
    monkeypatch.setattr(common, "REFERENCES_DIR", refs)
    monkeypatch.setattr(common, "SYNVOY_DIR", synvoy)
    monkeypatch.setattr(common, "DISCOVERY_DIRS", discovery_dirs)
    monkeypatch.setattr(merger, "DISCOVERY_DIRS", discovery_dirs)
    monkeypatch.setattr(panel, "PROJECT_ROOT", root)
    monkeypatch.setattr(panel, "METADATA_DIR", metadata)
    monkeypatch.setattr(summary, "PROJECT_ROOT", root)
    monkeypatch.setattr(summary, "METADATA_DIR", metadata)
    monkeypatch.setattr(figures, "RESULTS_DIR", results)
    monkeypatch.setattr(figures, "SYNVOY_DIR", synvoy)
    monkeypatch.setattr(phylo_input, "REFERENCES_DIR", refs)
    monkeypatch.setattr(phylo_input, "RESULTS_DIR", results)
    monkeypatch.setattr(tree_helper, "METADATA_DIR", metadata)
    monkeypatch.setattr(
        panel,
        "PANEL_SOURCES",
        (
            ("reference_innexin", metadata / "species_config.csv", "innexin"),
            ("reference_connexin", metadata / "species_config.csv", "connexin"),
            ("discovery_innexin", metadata / "not_annotated_innexin.txt", "innexin"),
        ),
    )

    return {
        "root": root,
        "metadata": metadata,
        "results": results,
        "refs": refs,
        "synvoy": synvoy,
        "discovery": discovery,
    }
