"""Tests for MMseqs region merge and Dmel truth loading."""

from __future__ import annotations

from pathlib import Path

import pipeline.innexin_locate as locate
from pipeline.dmel_innexin_benchmark import (
    GENE_TO_ACCESSION,
    intervals_overlap,
    load_truth_genes,
    score_recovery,
)
from pipeline.innexin_locate import LocatedLocus


def test_classify_ignores_terminal_stop():
    from pipeline.batch_innexin_curator_probe import classify_product

    # Complete ORF with terminal stop should still count as present
    prot = ("M" + "A" * 300 + "C" * 4 + "*")
    assert classify_product(0.5, prot) == "present"
    # Internal stop → not present
    bad = ("M" + "A" * 100 + "*" + "A" * 200 + "C" * 4 + "*")
    assert classify_product(0.5, bad) != "present"


def test_seed_regions_nms_keeps_nearby_paralogs():
    hits = [
        ("chrX", 1000, 3000, 1e-100),  # ogre-like
        ("chrX", 3500, 4500, 1e-40),  # Inx7-like (~500bp gap)
        ("chrX", 9000, 11000, 1e-200),  # Inx2-like
        ("chrX", 1200, 2800, 1e-50),  # redundant ogre hit
    ]
    regions = locate.seed_regions_nms(hits, suppress_bp=250, seed_pad_bp=100)
    assert len(regions) == 3
    starts = [r.start for r in regions]
    assert starts == sorted(starts)


def test_merge_hit_intervals_merges_nearby():
    hits = [
        ("chr1", 100, 200, 1e-20),
        ("chr1", 250, 400, 1e-10),  # within 50 kb gap
        ("chr1", 100_000, 100_100, 1e-5),
        ("chr2", 10, 20, 1e-8),
    ]
    regions = locate.merge_hit_intervals(hits, merge_gap_bp=50_000)
    assert len(regions) == 3
    assert regions[0].seqid == "chr1"
    assert regions[0].start == 100
    assert regions[0].end == 400
    assert regions[0].n_hits == 2
    assert regions[0].best_evalue == 1e-20


def test_intervals_overlap():
    assert intervals_overlap(1, 10, 8, 20)
    assert not intervals_overlap(1, 10, 20, 30)
    assert intervals_overlap(1, 10, 12, 20, pad=5)


def test_load_truth_genes():
    gff = Path(
        "project/data/annotations/innexin/Drosophila_melanogaster/"
        "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"
    )
    if not gff.exists():
        return
    truth = load_truth_genes(gff)
    assert len(truth) == 8
    assert {t.name for t in truth} == set(GENE_TO_ACCESSION)
    inx7 = next(t for t in truth if t.name == "Inx7")
    assert inx7.start == 6991612


def test_score_recovery_matches_overlap():
    truth = load_truth_genes() if Path(
        "project/data/annotations/innexin/Drosophila_melanogaster/"
        "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"
    ).exists() else []
    if not truth:
        return
    gene = next(t for t in truth if t.name == "Inx2")
    loci = [
        LocatedLocus(
            seqid=gene.seqid,
            start=gene.start + 10,
            end=gene.end - 10,
            strand="+",
            identity=0.9,
            target="Inx2",
            aa_length=367,
            cys=8,
            stops=0,
            verdict="present",
            exon_count=2,
            genomic_span=gene.end - gene.start,
            region_id="R0",
            method="test",
        )
    ]
    rows = score_recovery(truth, loci)
    by_gene = {r["gene"]: r for r in rows}
    assert by_gene["Inx2"]["recovered"] == "yes"
    assert by_gene["Inx3"]["recovered"] == "no"
