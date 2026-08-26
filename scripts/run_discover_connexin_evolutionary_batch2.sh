#!/usr/bin/env bash
# Connexin genome discovery for evolutionary clade species batch 2 (25 total).
set -uo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
export PATH="${HOME}/miniconda3/envs/synvoy_env/bin:${PATH}"
export PYTHONUNBUFFERED=1

mkdir -p project/data/genomes/connexin_evolutionary_batch2_discovery
mkdir -p project/results/connexin_evolutionary_batch2_discovery

python3 -u tools/discover_innexins.py \
  --family connexin \
  --species-txt project/metadata/not_annotated_connexin_evolutionary_batch2.txt \
  --genome-dir project/data/genomes/connexin_evolutionary_batch2_discovery \
  --results-dir project/results/connexin_evolutionary_batch2_discovery \
  --num-threads 4 \
  --chunk-workers 3 \
  --chunk-min-bp 1000000 \
  2>&1 | tee -a project/results/connexin_evolutionary_batch2_discovery/run.log
