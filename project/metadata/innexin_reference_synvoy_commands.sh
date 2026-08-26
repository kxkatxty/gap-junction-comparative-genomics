#!/usr/bin/env bash
# SynVoy runs: innexin reference panel only (species_config.csv).
# Insect and nematode clades are run separately — they do not share microsynteny.
set -euo pipefail
cd "$(dirname "$0")/../../SynVoy"
source "${HOME}/miniconda3/etc/profile.d/conda.sh" && conda activate synvoy_env

SYNVOY_ARGS=(--auto_params false --multi_profile false
  --n_flanking_genes 5 --adaptive_max_regions 3
  --enable_smith_waterman false --exon_level_search false
  -resume)

# Pro-mode paths for insect panel (skip NCBI re-download; genomes cached locally).
INSECT_STAGE="results/innexin_synvoy/refpanel_dmel_Inx2/staged_targets"
INSECT_PRO_ARGS=(--mode pro
  --home_genome results/innexin_synvoy/refpanel_dmel_Inx2/home_genome/home_genome/home_genome.fna
  --home_gff results/innexin_synvoy/refpanel_dmel_Inx2/home_genome/home_genome/home_genome.gff
  --home_species "Drosophila melanogaster"
  --target_genomes "${INSECT_STAGE}/*")

# Schistocerca excluded: large genome causes MMseqs OOM on this machine.
INSECT_TARGETS="Aedes aegypti,Anopheles gambiae"
INSECT_MAX_GENOMES=2

echo '=== insect panel: Inx2 (Q9V427) ==='
nextflow run main.nf -profile laptop_safe \
  --query results/innexin_synvoy/refpanel_dmel_Inx2/query/resolved_query/Q9V427.fasta \
  --outdir results/innexin_synvoy/refpanel_dmel_Inx2 \
  "${INSECT_PRO_ARGS[@]}" "${SYNVOY_ARGS[@]}"

echo '=== insect panel: shakB (P33085) ==='
nextflow run main.nf -profile laptop_safe \
  --query results/innexin_synvoy/refpanel_dmel_shakB/query/resolved_query/P33085.fasta \
  --outdir results/innexin_synvoy/refpanel_dmel_shakB \
  "${INSECT_PRO_ARGS[@]}" "${SYNVOY_ARGS[@]}"

echo '=== nematode panel: inx-2 (Q9U3K5) ==='
nextflow run main.nf -profile standard \
  --query_id Q9U3K5 \
  --target_species "Caenorhabditis briggsae" \
  --max_genomes 1 \
  --outdir results/innexin_synvoy/refpanel_cele_inx-2 \
  "${SYNVOY_ARGS[@]}"

echo 'Done. Plots in SynVoy/results/innexin_synvoy/refpanel_*/'
