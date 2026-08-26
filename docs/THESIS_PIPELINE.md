# Comparative genomics pipeline — gap-junction families

Goal: reproducible comparative annotation and evolution workflow for **innexins** and **connexins**.

## Module map

| Module | Status | Script |
|--------|--------|--------|
| species_panel_builder | **new** | `pipeline/species_panel_builder.py` |
| annotation_status_checker | exists | `tools/classify_family_annotation.py` |
| genome_quality_table_builder | exists | `tools/build_quality_table.py` |
| reference_sequence_downloader | exists | `tools/download_by_family.py`, `tools/download_gff3_from_species_list.py` |
| family_result_merger | **new** | `pipeline/family_result_merger.py` |
| gap_junction_discovery | exists | `tools/discover_innexins.py` (`--family innexin\|connexin`) |
| candidate_validator | **new** | `pipeline/candidate_validator.py` |
| exon_structure_extractor | exists | `tools/extract_exon_structures.py` |
| exon_map_plotter | exists | `tools/plot_exon_maps.py` |
| phylogeny_input_builder | **new** | `pipeline/phylogeny_input_builder.py` |
| phylogeny_runner | **new** | `pipeline/phylogeny_runner.py` |
| tree_annotation_helper | **new** | `pipeline/tree_annotation_helper.py` |
| summary_table_builder | **new** | `pipeline/summary_table_builder.py` |
| figure_export_pipeline | **new** | `pipeline/figure_export_pipeline.py` |
| case_study_synteny_runner | **new** | `pipeline/case_study_synteny_runner.py` |
| synteny_comparison_builder | **new** | `pipeline/synteny_comparison_builder.py` |

## Workflow phases

### Phase 1 — Dataset (run now)
```bash
./scripts/run_thesis_phase1.sh
```
Outputs:
- `project/metadata/species_panel_master.csv`
- `project/metadata/gap_junction_candidates_master.csv`
- `project/metadata/candidate_validation_summary.csv`
- `project/metadata/thesis_master_summary.csv`

### Phase 2 — Structure comparison
```bash
python3 tools/extract_exon_structures.py project/metadata/species_config.csv
python3 tools/plot_exon_maps.py
```

### Phase 3 — Discovery
```bash
python3 tools/discover_innexins.py --family connexin --species-txt <queue.txt> ...
python3 -m pipeline.family_result_merger --accepted-only
python3 -m pipeline.candidate_validator
```

### Phase 4 — Evolution
```bash
python3 -m pipeline.phylogeny_input_builder --family innexin --include-references
python3 -m pipeline.sequence_alignment_runner project/results/phylogeny/gap_junction_phylogeny_input.fasta
python3 -m pipeline.phylogeny_runner project/results/phylogeny/gap_junction_phylogeny_input.aln.fasta
python3 -m pipeline.tree_annotation_helper project/results/phylogeny/gap_junction_phylogeny_input.treefile
```

### Phase 5 — Synteny comparison (selected case studies)
```bash
# Compare existing SynVoy results
./scripts/run_thesis_synteny.sh --compare-only

# Run one gene, then compare
./scripts/run_thesis_synteny.sh --gene Inx5

# Figures export (includes synteny HTML plots)
python3 -m pipeline.figure_export_pipeline --clean
```
Outputs:
- `project/metadata/synteny_job_summary.csv`
- `project/metadata/synteny_genome_comparison.csv`
- `project/metadata/synteny_home_flanking_genes.csv`

## Key data locations

| Asset | Path |
|-------|------|
| Reference proteins | `project/data/references/` |
| Annotations (GFF) | `project/data/annotations/` |
| Innexin discovery | `project/results/innexin_discovery/` |
| Connexin discovery | `project/results/connexin_evolutionary_*_discovery/` |
| Exon structures | `project/results/exon_structures/` |
| SynVoy case studies | `SynVoy/results/innexin_synvoy/` |

## Current project state (snapshot)

- **Reference panel:** 31 species in `project/metadata/species_config.csv` (6 innexin + 25 connexin)
- **Discovery:** ~40 connexin species with accepted candidates; 5 innexin species with accepted hits
- **SynVoy:** 3 complete + 1 partial innexin case study (ogre, Inx2, shakB)
- **Exon maps:** generated for reference panel

## How we organised the work

1. **Clade-based discovery**, not random species lists
2. **SynVoy for case studies only**, not bulk screening
3. **One master summary table** for overview tables
4. **Modular scripts** — rerun filters without redoing discovery
