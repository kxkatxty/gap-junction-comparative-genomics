# Metadata

Species lists, annotation tables, and job queues. Master tables:

- `species_config.csv` — reference panel (innexin + connexin)
- `reference_species_list.txt` — plain species names
- `species_panel_master.csv` — combined panel
- `gap_junction_candidates_master.csv` — merged discovery candidates
- `thesis_master_summary.csv` — overview numbers (filename kept for compatibility)

## Groups (by filename prefix)

| Prefix / pattern | Meaning |
|------------------|---------|
| `species_*` | Panel, catalogs, discovery queues |
| `not_annotated_*` | Species lacking family annotation |
| `family_annotation_*` | Classification status/summaries |
| `target_clades_*` | Clade targets for discovery |
| `*_synvoy_queue*` | SynVoy job queues |
| `synteny_*` | Synteny comparison tables |
| `connexin_*_parts`, `*_batches` | Parallel-job split files |
