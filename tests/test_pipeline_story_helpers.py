from __future__ import annotations

import csv
from pathlib import Path

import pipeline.gff_synteny_neighbor_summary as neigh
import pipeline.new_species_gallery_builder as gallery
import pipeline.innexin_subfamily_hypothesis as hypo
import pipeline.merge_curator_probe_into_master as merge
import pipeline.batch_innexin_curator_probe as curator
import pipeline.innexin_reference_panel_builder as refpanel
import pipeline.synteny_case_study_1_dmel_cluster as case1
import pipeline.phylogenetic_story_builder as story
import pipeline.gap_junction_path_builder as path
import pipeline.synteny_phylo_guided_inx_cluster as phylo_guided


class TestPhyloGuided:
    def test_phylo_panels_returns_list(self):
        panels = phylo_guided.phylo_panels()
        assert isinstance(panels, list)
        for item in panels:
            assert len(item) == 3
            assert isinstance(item[0], str)
            assert isinstance(item[2], list)


class TestNeighborSummary:
    def test_named_gene_and_panel_meta(self):
        assert neigh._is_named_gene("GJA1") is True
        assert neigh._is_named_gene("LOC123") is False
        assert neigh._is_named_gene("MIR-1") is False
        assert neigh._panel_meta("innexin_Inx2_cross_species") == ("innexin", "Inx2")
        assert neigh._panel_meta("connexin_GJA1_cross_species") == ("connexin", "GJA1")
        assert neigh._panel_meta("other") == ("other", "other")

    def test_conservation_tier(self):
        assert neigh._conservation_tier(3, 0, 5) == "high"
        assert neigh._conservation_tier(1, 0, 5) == "moderate"
        assert neigh._conservation_tier(0, 3, 5) == "moderate"
        assert neigh._conservation_tier(0, 1, 5) == "low"
        assert neigh._conservation_tier(0, 0, 5) == "none"

    def test_analyze_manifest(self, tmp_path: Path):
        manifest = tmp_path / "manifest.csv"
        with manifest.open("w", encoding="utf-8", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=["panel", "species", "name", "is_goi"],
            )
            writer.writeheader()
            writer.writerows(
                [
                    {"panel": "innexin_Inx2_cross_species", "species": "Dmel", "name": "Inx2", "is_goi": "yes"},
                    {"panel": "innexin_Inx2_cross_species", "species": "Dmel", "name": "ogre", "is_goi": "no"},
                    {"panel": "innexin_Inx2_cross_species", "species": "Agam", "name": "ogre", "is_goi": "no"},
                    {"panel": "innexin_Inx2_cross_species", "species": "Agam", "name": "LOC99", "is_goi": "no"},
                ]
            )
        detail, summary, conclusions = neigh.analyze_manifest(manifest)
        assert any(r["neighbor_gene"] == "ogre" for r in detail)
        assert summary
        assert isinstance(conclusions, dict)


class TestGalleryHelpers:
    def test_is_new_and_display(self):
        assert gallery._is_new({"prior_discovery": "not_run"}) is True
        assert gallery._is_new({"prior_discovery": "accepted=2", "prior_accepted": "2"}) is False
        assert gallery._is_new({"prior_discovery": "accepted=0", "prior_accepted": "0"}) is True
        assert gallery._display_name("Adineta_vaga") == "Adineta vaga"

    def test_bar_and_build_html(self):
        assert "width:50%" in gallery._bar(0.5)
        species = [
            {
                "name": "Adineta vaga",
                "clade": "Rotifera",
                "habit": "bdelloid",
                "status": "new",
                "status_label": "new",
                "is_new": True,
                "present": 2,
                "best_aa": 350,
                "best_id": 0.42,
                "best_hit": "unc-9",
                "best_locus": "cand_1",
                "prior": "not_run",
                "note": "",
            }
        ]
        html = gallery.build_html(species)
        assert "Adineta vaga" in html
        assert "Rotifera" in html


class TestSubfamilyHypothesis:
    def test_family_stem(self):
        assert hypo._family_stem("ABC1") == "ABC"
        assert hypo._family_stem("CG12345") == "CG12345"
        assert hypo._family_stem("lncRNA-x") == ""
        assert hypo._family_stem("mir-1") == ""

    def test_write_hypothesis_tables(self, tmp_path: Path):
        hypo.write_hypothesis_tables(tmp_path)
        assert (tmp_path / "subfamily_definitions.csv").exists()
        assert (tmp_path / "tip_subfamily_assignments.csv").exists()
        defs = list(csv.DictReader((tmp_path / "subfamily_definitions.csv").open(encoding="utf-8")))
        assert len(defs) >= 4


class TestMergeCurator:
    def test_organism(self):
        assert merge._organism("Adineta_vaga") == "Adineta vaga"

    def test_curator_rows_ranking(self, tmp_path: Path, monkeypatch):
        summary = tmp_path / "summary.csv"
        summary.write_text(
            "species,status,prior_discovery,prior_accepted\n"
            "New_sp,searched,not_run,\n"
            "Old_sp,searched,accepted=2,2\n",
            encoding="utf-8",
        )
        hits_new = tmp_path / "New_sp"
        hits_new.mkdir()
        (hits_new / "hits.csv").write_text(
            "verdict,aa_length,seqid,start,end,strand,identity,target,cys,stops\n"
            "present,350,chr1,1,1000,+,0.4,unc-9,6,0\n"
            "fragmentary,200,chr1,2000,2500,+,0.3,unc-7,4,0\n"
            "weak,120,chr2,1,400,+,0.25,inx-2,2,0\n",
            encoding="utf-8",
        )
        monkeypatch.setattr(merge, "CURATOR_SUMMARY", summary)
        monkeypatch.setattr(merge, "CURATOR_ROOT", tmp_path)
        rows = merge.curator_rows(include_known=False)
        assert len(rows) == 2
        assert rows[0]["rank_category"] == "curator_present_innexin"
        assert rows[1]["rank_category"] == "weak_manual_review"


class TestCuratorProbe:
    def test_classify_product(self):
        long = "C" * 4 + "A" * 300
        assert curator.classify_product(0.25, long) == "present"
        assert curator.classify_product(0.25, "C" * 4 + "A" * 200) == "fragmentary"
        assert curator.classify_product(0.25, "A" * 120) == "weak"
        assert curator.classify_product(0.1, "A" * 400) == "noise"

    def test_find_genome_fasta(self, tmp_path: Path):
        (tmp_path / "small.fna").write_text("A" * 10, encoding="utf-8")
        (tmp_path / "big.decompressed.fna").write_text("A" * 100, encoding="utf-8")
        chosen = curator.find_genome_fasta(tmp_path)
        assert chosen is not None
        assert "decompressed" in chosen.name

    def test_prior_discovery(self, tmp_path: Path, monkeypatch):
        monkeypatch.setattr(curator, "DISCOVERY_ROOT", tmp_path)
        assert curator.prior_discovery("Missing") == ("not_run", "")
        sp = tmp_path / "Sp"
        sp.mkdir()
        (sp / "discovery_summary.json").write_text('{"accepted_count": 2, "candidate_count": 5}', encoding="utf-8")
        status, accepted = curator.prior_discovery("Sp")
        assert "accepted=2" in status
        assert accepted == "2"

    def test_parse_miniprot_models(self, tmp_path: Path):
        gff = tmp_path / "out.gff"
        gff.write_text(
            "##gff-version 3\n"
            "chr1\tminiprot\tmRNA\t10\t100\t.\t+\t.\tID=m1;Identity=0.4;Target=Inx2\n"
            "chr1\tminiprot\tCDS\t10\t40\t.\t+\t0\tParent=m1\n"
            "chr1\tminiprot\tCDS\t50\t100\t.\t+\t0\tParent=m1\n",
            encoding="utf-8",
        )
        models = curator.parse_miniprot_models(gff)
        assert len(models) == 1
        assert models[0]["identity"] == 0.4
        assert len(models[0]["cds"]) == 2


class TestReferencePanel:
    def test_load_and_target_species(self, tmp_path: Path):
        cfg = tmp_path / "species_config.csv"
        cfg.write_text(
            "family,organism,reviewed_only\n"
            "innexin,Drosophila melanogaster,true\n"
            "innexin,Caenorhabditis elegans,true\n"
            "connexin,Homo sapiens,true\n",
            encoding="utf-8",
        )
        species = refpanel.load_innexin_species(cfg)
        assert species == ["Caenorhabditis elegans", "Drosophila melanogaster"]
        targets = refpanel.target_species("Drosophila melanogaster", species)
        assert "Caenorhabditis elegans" in targets
        assert "Drosophila melanogaster" not in targets

    def test_build_synvoy_script(self):
        panel = ["Drosophila melanogaster", "Aedes aegypti", "Caenorhabditis elegans"]
        script = refpanel.build_synvoy_script(panel)
        assert "set -euo pipefail" in script
        assert "synvoy_env" in script

    def test_panel_rows_with_sandbox(self, pipeline_sandbox):
        rows = refpanel.panel_rows(["Drosophila melanogaster"])
        assert rows[0]["organism"] == "Drosophila melanogaster"
        assert rows[0]["clade"] == "insect"
        assert int(rows[0]["reference_protein_count"]) >= 1


class TestCase1AndStoryAndPath:
    def test_bed_and_plot_paths(self):
        bed = case1.bed_path("dmel_inx2", "synteny_block_1")
        plot = case1.plot_path("dmel_inx2", "synteny_block_1")
        assert bed.name.endswith(".bed")
        assert "plot_inputs_synteny_block_1" in str(bed)
        assert plot.name.endswith("_synteny_plot.html")

    def test_phylogenetic_story_build(self, tmp_path: Path, monkeypatch):
        src_dir = tmp_path / "src"
        src_dir.mkdir()
        figures = []
        for name in ("a.png", "b.png"):
            p = src_dir / name
            p.write_bytes(b"\x89PNG\r\n\x1a\n")
            figures.append((name, p))
        monkeypatch.setattr(story, "STORY_FIGURES", figures)
        monkeypatch.setattr(story, "RESULTS_DIR", tmp_path)
        out = tmp_path / "story"
        story.build(out)
        assert (out / "figures" / "a.png").exists()
        assert (out / "figures" / "b.png").exists()
        html = (out / "index.html").read_text(encoding="utf-8")
        assert "Pannexin bridge" in html
        assert "How gap junctions diverged" in html

    def test_story_figure_manifest_includes_pannexin(self):
        names = [n for n, _ in story.STORY_FIGURES]
        assert "05_pannexin_bridge.png" in names
        assert "06_panx_vs_inx.png" in names
        assert "07c_panx_iqtree.png" in names
        assert "08_connexin_synteny.png" in names

    def test_gap_junction_path_plots(self, tmp_path: Path):
        for fn in (path.plot_two_cells, path.plot_topology, path.plot_hexamers, path.plot_kingdoms):
            out = tmp_path / f"{fn.__name__}.png"
            fn(out)
            assert out.exists()
            assert out.stat().st_size > 0
