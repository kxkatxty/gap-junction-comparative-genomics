#!/usr/bin/env bash
# FOXP2-style microsynteny gene-order plots (NCBI GFF, no SynVoy).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
export MPLCONFIGDIR="${ROOT}/project/results/synteny_foxp2_style/.mpl_cache"
mkdir -p "$MPLCONFIGDIR"
python3 -m pipeline.foxp2_style_synteny --outdir project/results/synteny_foxp2_style
echo "Gallery: project/results/synteny_foxp2_style/index.html"
