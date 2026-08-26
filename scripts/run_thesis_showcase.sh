#!/usr/bin/env bash
# Build the HTML showcase gallery from pipeline figures.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
export MPLCONFIGDIR="${ROOT}/project/results/showcase/.mpl_cache"
mkdir -p "$MPLCONFIGDIR"

echo "=== Step 1: Exon & protein feature maps (existing pipeline) ==="
python3 tools/extract_exon_structures.py 2>/dev/null || true
python3 tools/plot_exon_maps.py 2>/dev/null || true

echo "=== Step 2: Scientific showcase (35+ comparison plots) ==="
python3 -m pipeline.scientific_showcase_builder --outdir project/results/showcase

echo "=== Step 3: Copy key existing plots into showcase ==="
DEST="project/results/showcase/imported"
mkdir -p "$DEST/exon_maps" "$DEST/protein_features" "$DEST/synteny_case_studies"
for f in project/results/exon_structures/plots/*.png; do
  [ -f "$f" ] && cp -f "$f" "$DEST/exon_maps/"
done
for f in project/results/exon_structures/plots/protein_feature_maps/*.png; do
  [ -f "$f" ] && cp -f "$f" "$DEST/protein_features/"
done
for f in project/results/synteny_case_studies/case*/figures/*.html; do
  [ -f "$f" ] && cp -f "$f" "$DEST/synteny_case_studies/" 2>/dev/null || true
done

echo "=== Done ==="
echo "Gallery: project/results/showcase/index.html"
echo "Plots:   project/results/showcase/plot_manifest.csv"
