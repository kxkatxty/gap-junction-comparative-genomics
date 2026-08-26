from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import collect_species_by_clade as csc


class TestCladeIO:
    def test_read_clade_targets(self, tmp_path: Path):
        path = tmp_path / "clades.csv"
        path.write_text(
            "clade,ncbi_search_term,notes,priority\n"
            "Rotifera,Rotifera[Organism],rotifers,4\n"
            "Annelida,Annelida[Organism],worms,5\n",
            encoding="utf-8",
        )
        targets = csc.read_clade_targets(path)
        assert [t.clade for t in targets] == ["Rotifera", "Annelida"]
        assert targets[0].ncbi_search_term == "Rotifera[Organism]"


class TestAssemblyHelpers:
    def test_clade_assembly_search_term_adds_latest_filter(self):
        assert "latest[filter]" in csc.clade_assembly_search_term("Rotifera[Organism]")

    def test_best_assembly_per_organism_picks_ranked(self):
        records = [
            {
                "Organism": "Artemia franciscana",
                "AssemblyAccession": "GCA_1",
                "FtpPath_RefSeq": "",
                "RefSeq_category": "",
                "AssemblyStatus": "Contig",
            },
            {
                "Organism": "Artemia franciscana",
                "AssemblyAccession": "GCF_2",
                "FtpPath_RefSeq": "ftp://example",
                "RefSeq_category": "reference genome",
                "AssemblyStatus": "Chromosome",
            },
        ]
        best = csc.best_assembly_per_organism(records)
        assert best["Artemia franciscana"]["AssemblyAccession"] == "GCF_2"


class TestPriority:
    def test_discovery_priority_for_unannotated_species(self):
        candidate = csc.SpeciesCandidate(
            clade="Rotifera",
            organism="Rotaria sordida",
            assembly_accession="GCF_123",
            assembly_level="Scaffold",
            has_gff_annotation="yes",
            innexin_annotation_status="no_related_entry_found",
            already_checked="yes",
        )
        scored = csc.compute_priority(candidate)
        assert scored.recommended_pipeline == "discovery"
        assert scored.priority_tier in {"high", "medium"}

    def test_annotated_species_routes_to_annotated_pipeline(self):
        candidate = csc.SpeciesCandidate(
            clade="Crustacea",
            organism="Daphnia pulex",
            assembly_accession="GCF_123",
            assembly_level="Chromosome",
            innexin_annotation_status="annotated_family_member_found",
        )
        scored = csc.compute_priority(candidate)
        assert scored.recommended_pipeline == "annotated"


class TestFiltering:
    def test_filter_candidates_exclude_model_and_checked(self):
        candidates = [
            csc.SpeciesCandidate(
                clade="Rotifera",
                organism="Rotaria macrura",
                assembly_accession="GCF_1",
                assembly_level="Scaffold",
                already_checked="yes",
                model_species="no",
            ),
            csc.SpeciesCandidate(
                clade="Rotifera",
                organism="Drosophila melanogaster",
                assembly_accession="GCF_2",
                assembly_level="Chromosome",
                already_checked="no",
                model_species="yes",
            ),
        ]
        kept = csc.filter_candidates(
            candidates,
            exclude_model_species=True,
            exclude_already_checked=True,
            min_assembly_level="scaffold",
            require_gff=False,
            pipeline=None,
        )
        assert kept == []

    def test_clean_organism_name(self):
        assert csc.clean_organism_name("Trichoplax adhaerens (placozoans)") == "Trichoplax adhaerens"

    def test_apply_per_clade_limit(self):
        candidates = [
            csc.SpeciesCandidate(clade="Rotifera", organism="A", priority_score=10),
            csc.SpeciesCandidate(clade="Rotifera", organism="B", priority_score=20),
            csc.SpeciesCandidate(clade="Annelida", organism="C", priority_score=30),
        ]
        limited = csc.apply_per_clade_limit(candidates, 1)
        assert {item.organism for item in limited} == {"B", "C"}


class TestCollection:
    def test_collect_candidates_for_clade(self):
        target = csc.CladeTarget("Rotifera", "Rotifera[Organism]")
        records = [
            {
                "Organism": "Rotaria macrura",
                "AssemblyAccession": "GCF_1",
                "AssemblyStatus": "Scaffold",
                "FtpPath_RefSeq": "ftp://example",
            }
        ]
        with patch.object(csc, "fetch_assemblies_for_clade", return_value=records):
            with patch.object(csc, "has_gff_for_assembly", return_value=False):
                candidates = csc.collect_candidates_for_clade(
                    target,
                    master={},
                    summary={},
                    model_species=set(),
                    check_gff=False,
                    retmax=10,
                    request_delay_sec=0,
                )
        assert len(candidates) == 1
        assert candidates[0].organism == "Rotaria macrura"
        assert candidates[0].clade == "Rotifera"
