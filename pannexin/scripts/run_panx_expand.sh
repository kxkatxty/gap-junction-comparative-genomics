#!/usr/bin/env bash
# Light pannexin rebuild only — no SynVoy, no genome miniprot, no tree recompute.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"

export PATH="${HOME}/miniconda3/envs/synvoy_env/bin:${PATH}"
export PANX_THREADS="${PANX_THREADS:-1}"
export PANX_HEAVY_PHYLO=0
export PANX_RUN_MINIPROT=0

echo "Light pannexin rebuild (1 thread; SynVoy deferred; no miniprot)"

python3 -m pipeline.pannexin_panel_builder
# UniProt/NCBI status for panel species (network). Reuse with: PANX_REUSE_ANNOTATION=1
if [[ "${PANX_REUSE_ANNOTATION:-0}" == "1" ]]; then
  python3 -m pipeline.pannexin_annotation_status_builder --reuse
else
  python3 -m pipeline.pannexin_annotation_status_builder
fi
python3 -m pipeline.pannexin_clade_comparison_builder
python3 -m pipeline.pannexin_similarity_builder
python3 -m pipeline.panx_vs_inx_builder
python3 -m pipeline.pannexin_phylogeny_builder   # reuses existing tree; renders figure
# Exon notes need network once; reuse with PANX_REUSE_EXONS=1
if [[ "${PANX_REUSE_EXONS:-0}" == "1" ]]; then
  python3 -m pipeline.pannexin_exon_notes_builder --reuse
else
  python3 -m pipeline.pannexin_exon_notes_builder
fi
python3 -m pipeline.pannexin_curator_builder     # status table only
python3 -m pipeline.pannexin_path_builder
python3 -m pipeline.pannexin_insights_builder

echo
echo "Open:"
echo "  project/results/pannexin_insights/index.html"
echo "  project/results/pannexin_phylogeny/index.html"
echo "  project/results/pannexin_path/index.html"
echo
echo "SynVoy / genome probes: not part of this script — do later, separately."
