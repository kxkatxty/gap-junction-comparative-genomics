#!/usr/bin/env bash
# Compare the 6 innexin reference species: summary table + reference-only phylogeny.
set -euo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"

echo "=== innexin reference panel table ==="
python3 -m pipeline.innexin_reference_panel_builder

echo "=== reference-only phylogeny FASTA ==="
python3 -m pipeline.phylogeny_input_builder \
  --family innexin \
  --include-references \
  --output project/results/phylogeny/innexin_reference_panel.fasta

if command -v mafft >/dev/null 2>&1; then
  python3 -m pipeline.sequence_alignment_runner \
    project/results/phylogeny/innexin_reference_panel.fasta \
    --output project/results/phylogeny/innexin_reference_panel.aln.fasta
fi

if command -v iqtree >/dev/null 2>&1 || command -v iqtree2 >/dev/null 2>&1; then
  python3 -m pipeline.phylogeny_runner \
    project/results/phylogeny/innexin_reference_panel.aln.fasta \
    --prefix project/results/phylogeny/innexin_reference_panel
fi

echo ""
echo "Outputs:"
echo "  project/metadata/innexin_reference_panel.csv"
echo "  project/metadata/innexin_reference_synvoy_commands.sh  (optional SynVoy)"
echo "  project/results/phylogeny/innexin_reference_panel.*"
echo ""
echo "To run SynVoy across reference insects/nematodes separately:"
echo "  ./project/metadata/innexin_reference_synvoy_commands.sh"
