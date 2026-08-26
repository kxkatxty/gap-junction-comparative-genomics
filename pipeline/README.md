# Pipeline modules

Run as `python3 -m pipeline.<module>` from the repo root.

## HTML builders

- `gap_junction_path_builder` — main path
- `phylogenetic_story_builder`, `new_species_gallery_builder`
- `innexin_*` / `connexin_*` clade, similarity, insights
- `inx_vs_cnx_builder` — family contrast
- `family_comparison_builder` — side-by-side metrics
- `scientific_showcase_builder`

## Core steps

- `species_panel_builder`, `family_result_merger`, `candidate_validator`
- `summary_table_builder`, `figure_export_pipeline`
- `phylogeny_input_builder`, `sequence_alignment_runner`, `phylogeny_runner`

## Synteny

- `gff_microsynteny`, `foxp2_style_synteny`
- `case_study_synteny_runner`, `synteny_comparison_builder`
- `synteny_case_study_1_*` … `3_*`, `synteny_phylo_guided_inx_cluster`

Shared paths: `pipeline/common.py`
