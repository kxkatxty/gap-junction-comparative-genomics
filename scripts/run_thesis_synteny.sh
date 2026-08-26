#!/usr/bin/env bash
# Synteny case studies: run selected SynVoy jobs + build comparison tables.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"

GENES=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --gene) GENES+=("$2"); shift 2 ;;
    --compare-only) COMPARE_ONLY=1; shift ;;
    *) echo "Usage: $0 [--compare-only] [--gene GENE ...]"; exit 1 ;;
  esac
done

if [[ "${COMPARE_ONLY:-0}" != "1" ]]; then
  if [[ ${#GENES[@]} -eq 0 ]]; then
    echo "No --gene given; building comparison tables from existing SynVoy results only."
  else
    ARGS=()
    for g in "${GENES[@]}"; do ARGS+=(--gene "$g"); done
    python3 -m pipeline.case_study_synteny_runner "${ARGS[@]}"
  fi
fi

python3 -m pipeline.synteny_comparison_builder
python3 -m pipeline.figure_export_pipeline --clean

echo "Done. See project/metadata/synteny_*.csv and project/results/figures/thesis/"
