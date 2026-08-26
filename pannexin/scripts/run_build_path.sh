#!/usr/bin/env bash
# Build the guided pannexin reading path (HTML + starter figures).
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
python3 -m pipeline.pannexin_path_builder
echo "Open: project/results/pannexin_path/index.html"
