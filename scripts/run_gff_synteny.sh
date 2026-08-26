#!/usr/bin/env bash
# GFF-based microsynteny for gap-junction genes (no SynVoy).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
export MPLCONFIGDIR="${ROOT}/project/results/synteny_gff/.mpl_cache"
mkdir -p "$MPLCONFIGDIR"

MODE="${1:-all}"
python3 -m pipeline.gff_microsynteny --outdir project/results/synteny_gff --mode "$MODE"
echo "Gallery: project/results/synteny_gff/index.html"
