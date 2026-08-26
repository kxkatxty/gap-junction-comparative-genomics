#!/usr/bin/env bash
# Phase 2–4 helpers: exon maps, phylogeny input, optional alignment/tree.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"

FAMILY="${1:-innexin}"

echo "=== exon structure plots (if not yet generated) ==="
python3 tools/plot_exon_maps.py || true

echo "=== phylogeny input FASTA ==="
python3 -m pipeline.phylogeny_input_builder --family "$FAMILY" --include-references

if command -v mafft >/dev/null 2>&1; then
  echo "=== MAFFT alignment ==="
  python3 -m pipeline.sequence_alignment_runner "project/results/phylogeny/gap_junction_phylogeny_input.fasta"
else
  echo "SKIP alignment: mafft not installed"
fi

if command -v iqtree2 >/dev/null 2>&1 || command -v iqtree >/dev/null 2>&1; then
  echo "=== IQ-TREE ==="
  python3 -m pipeline.phylogeny_runner "project/results/phylogeny/gap_junction_phylogeny_input.aln.fasta"
  python3 -m pipeline.tree_annotation_helper "project/results/phylogeny/gap_junction_phylogeny_input.aln.treefile"
else
  echo "SKIP tree: iqtree not installed"
fi

echo "Done."
