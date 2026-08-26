#!/usr/bin/env bash
# Retry connexin discovery for batch2 species with failed genome downloads.
set -uo pipefail

source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
export PATH="${HOME}/miniconda3/envs/synvoy_env/bin:${PATH}"
export PYTHONUNBUFFERED=1

mkdir -p project/data/genomes/connexin_evolutionary_batch2_discovery
mkdir -p project/results/connexin_evolutionary_batch2_discovery

python3 -u tools/discover_innexins.py \
  --family connexin \
  --species-txt project/metadata/connexin_evolutionary_batch2_retry.txt \
  --genome-dir project/data/genomes/connexin_evolutionary_batch2_discovery \
  --results-dir project/results/connexin_evolutionary_batch2_discovery \
  --force-redownload \
  --num-threads 4 \
  --chunk-workers 3 \
  --chunk-min-bp 1000000 \
  2>&1 | tee -a project/results/connexin_evolutionary_batch2_discovery/retry.log
