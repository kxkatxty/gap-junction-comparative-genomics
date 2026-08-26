#!/usr/bin/env bash
# Phase 1: stabilize dataset — panel, candidates, validation, master summary.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"

echo "=== species panel ==="
python3 -m pipeline.species_panel_builder

echo "=== merge discovery candidates ==="
python3 -m pipeline.family_result_merger

echo "=== candidate validation summary ==="
python3 -m pipeline.candidate_validator

echo "=== master thesis summary ==="
python3 -m pipeline.summary_table_builder

echo "Done. See project/metadata/thesis_master_summary.csv"
