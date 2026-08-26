#!/usr/bin/env bash
# Connexin genome discovery for evolutionary clade species (25 total).
set -uo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
export PATH="${HOME}/miniconda3/envs/synvoy_env/bin:${PATH}"

python3 tools/discover_innexins.py \
  --family connexin \
  --species-txt project/metadata/not_annotated_connexin_evolutionary.txt \
  --genome-dir project/data/genomes/connexin_evolutionary_discovery \
  --results-dir project/results/connexin_evolutionary_discovery \
  2>&1 | tee project/results/connexin_evolutionary_discovery/run.log
