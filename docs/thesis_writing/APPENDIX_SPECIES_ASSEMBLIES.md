# Appendix: species, assemblies, dataset roles, and counts

Auto-generated from live project metadata (2026-09-09).

**Primary machine-readable file:** [`APPENDIX_SPECIES_ASSEMBLIES.csv`](APPENDIX_SPECIES_ASSEMBLIES.csv)

Also copied to `project/results/thesis_appendix/appendix_species_assemblies.csv`.

## Column definitions

| Column | Meaning |
|--------|---------|
| `family` | innexin / connexin / pannexin |
| `organism` | Species name |
| `assembly_accession` | NCBI GCA_/GCF_ when recoverable; may be blank for UniProt-only protein panels |
| `assembly_level` | NCBI assembly level when known |
| `dataset_roles` | Semicolon-separated roles in this thesis |
| `n_sequences_analytical_clade_table` | Sequences for that species in the main analytical table |
| `n_curator_present_floor` | Innexin curator present-locus floor (not a census) |

## Panel totals

- Rows: **112** (innexin 57, connexin 38, pannexin 17)
- Innexin analytical sequence sum: **325** (expect 325)
- Connexin analytical sum: **127** (expect 127)
- Pannexin panel sum: **49** (expect 49)
- Assembly accession filled: innexin 57/57, connexin 38/38, pannexin 9/17

## Innexin true-rescue species

| organism | assembly_accession | n_curator_present_floor | prior_discovery |
|----------|-------------------:|------------------------:|-----------------|
| Abra alba | GCA_979025805.1 | 11 | accepted=0;candidates=1 |
| Brachionus manjavacas | GCA_019054805.1 | 5 | accepted=0;candidates=1 |
| Rotaria macrura | GCA_900239685.1 | 11 | accepted=0;candidates=3 |

## Limits

- Accessions are best-effort merges across genome paths, GFF download report, quality table, and clade catalogs.
- Pannexin rows are UniProt protein-based; genome assembly is often not applicable (8/17 lack accession).
- Curator present counts are floors, not paralog censuses.
- GFF annotation *release date* is not always stored separately from the assembly accession.

