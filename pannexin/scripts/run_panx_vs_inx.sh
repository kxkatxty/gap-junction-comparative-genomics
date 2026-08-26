#!/usr/bin/env bash
# Build pannexin panel, within-family similarity, and panx×inx contrast.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"

python3 -m pipeline.pannexin_panel_builder
python3 -m pipeline.pannexin_similarity_builder
python3 -m pipeline.panx_vs_inx_builder
python3 -m pipeline.pannexin_path_builder

echo
echo "Open:"
echo "  project/results/panx_vs_inx/index.html"
echo "  project/results/pannexin_similarity/index.html"
echo "  project/results/pannexin_path/index.html"
