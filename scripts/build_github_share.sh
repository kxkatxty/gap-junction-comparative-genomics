#!/usr/bin/env bash
# Build a lean folder suitable for GitHub (~tens of MB, not hundreds of GB).
# Copies websites, builders, tools, metadata, and light curator summaries only.
set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
OUT="${1:-$ROOT/github_share}"

echo "Building github_share → $OUT"
# Preserve an existing git history (do not wipe .git).
if [[ -d "$OUT/.git" ]]; then
  TMP_GIT=$(mktemp -d)
  mv "$OUT/.git" "$TMP_GIT/"
  rm -rf "$OUT"
  mkdir -p "$OUT"
  mv "$TMP_GIT/.git" "$OUT/"
  rmdir "$TMP_GIT"
else
  rm -rf "$OUT"
  mkdir -p "$OUT"
fi

RSYNC=(rsync -a --delete
  --exclude '.DS_Store'
  --exclude '__pycache__/'
  --exclude '.pytest_cache/'
  --exclude '*.pyc'
  --exclude 'work/'
  --exclude 'tmp/'
  --exclude '_work/'
  --exclude '_cache/'
  --exclude 'mmseqs_tmp/'
  --exclude '*.log'
  --exclude '*.gz'
  --exclude '*.ckp'
)

# --- top-level docs / config ---
cp "$ROOT/README.md" "$OUT/README_FULL_REPO.md"
cp "$ROOT/pytest.ini" "$OUT/"
cp "$ROOT/requirements-dev.txt" "$OUT/"
cp "$ROOT/.gitignore" "$OUT/"

# --- code ---
"${RSYNC[@]}" "$ROOT/pipeline/" "$OUT/pipeline/"
"${RSYNC[@]}" "$ROOT/tools/" "$OUT/tools/"
"${RSYNC[@]}" "$ROOT/scripts/" "$OUT/scripts/"
"${RSYNC[@]}" "$ROOT/tests/" "$OUT/tests/"
"${RSYNC[@]}" "$ROOT/docs/" "$OUT/docs/"

# --- metadata (species panels, summary tables) ---
mkdir -p "$OUT/project"
"${RSYNC[@]}" "$ROOT/project/metadata/" "$OUT/project/metadata/"

# --- HTML result sites (no discovery dumps, no genomes) ---
RESULTS_KEEP=(
  gap_junction_path
  phylogenetic_story
  inx_vs_cnx
  family_comparison
  innexin_insights
  connexin_insights
  innexin_clade_comparison
  connexin_clade_comparison
  innexin_similarity
  connexin_similarity
  innexin_subfamilies
  new_species_gallery
  sources
  showcase
  synteny_gff
  synteny_case_studies
  synteny_phylo_guided
  synteny_foxp2_style
  exon_structures
  figures
  phylogeny
)

mkdir -p "$OUT/project/results"
if [[ -f "$ROOT/project/results/README.md" ]]; then
  cp "$ROOT/project/results/README.md" "$OUT/project/results/README.md"
fi

for name in "${RESULTS_KEEP[@]}"; do
  src="$ROOT/project/results/$name"
  if [[ -e "$src" ]]; then
    echo "  + results/$name"
    "${RSYNC[@]}" "$src" "$OUT/project/results/"
  fi
done

# Curator reproducibility pack (small): summaries + per-species hits/products + query pack
mkdir -p "$OUT/project/results/gene_curator_probe"
shopt -s nullglob
for f in "$ROOT/project/results/gene_curator_probe"/*.{md,csv,tsv,txt,fasta}; do
  cp "$f" "$OUT/project/results/gene_curator_probe/"
done
# Per-species: hits.csv + product FASTAs only (skip locate work dirs / large GFF dumps)
for sp_dir in "$ROOT/project/results/gene_curator_probe"/*/; do
  [[ -d "$sp_dir" ]] || continue
  sp=$(basename "$sp_dir")
  # skip non-species dirs if any
  [[ -f "$sp_dir/hits.csv" ]] || continue
  dest="$OUT/project/results/gene_curator_probe/$sp"
  mkdir -p "$dest"
  cp "$sp_dir/hits.csv" "$dest/"
  [[ -f "$sp_dir/present_products.fasta" ]] && cp "$sp_dir/present_products.fasta" "$dest/"
  [[ -f "$sp_dir/genome_direct_products.fasta" ]] && cp "$sp_dir/genome_direct_products.fasta" "$dest/"
  # light miniprot GFF if small (<2 MB)
  if [[ -f "$sp_dir/miniprot.gff" ]]; then
    sz=$(wc -c < "$sp_dir/miniprot.gff" | tr -d ' ')
    if [[ "$sz" -lt 2000000 ]]; then
      cp "$sp_dir/miniprot.gff" "$dest/"
    fi
  fi
done
shopt -u nullglob
echo "  + results/gene_curator_probe (summaries + hits/products)"

# --- reference protein panels (small; required for reproducibility) ---
mkdir -p "$OUT/project/data/references"
for fam in innexins connexins; do
  if [[ -d "$ROOT/project/data/references/$fam" ]]; then
    echo "  + data/references/$fam"
    "${RSYNC[@]}" "$ROOT/project/data/references/$fam" "$OUT/project/data/references/"
  fi
done
# pannexin refs live under pannexin/ (copied below)

# --- pannexin workspace (sites + scripts; empty data placeholders only) ---
mkdir -p "$OUT/pannexin/project/results" "$OUT/pannexin/project/metadata" \
  "$OUT/pannexin/project/data/references/pannexins" \
  "$OUT/pannexin/project/data/annotations/pannexin"
"${RSYNC[@]}" "$ROOT/pannexin/pipeline/" "$OUT/pannexin/pipeline/"
"${RSYNC[@]}" "$ROOT/pannexin/scripts/" "$OUT/pannexin/scripts/"
"${RSYNC[@]}" "$ROOT/pannexin/tests/" "$OUT/pannexin/tests/" 2>/dev/null || true
"${RSYNC[@]}" "$ROOT/pannexin/docs/" "$OUT/pannexin/docs/" 2>/dev/null || true
cp "$ROOT/pannexin/README.md" "$OUT/pannexin/" 2>/dev/null || true
cp "$ROOT/pannexin/pytest.ini" "$OUT/pannexin/" 2>/dev/null || true
cp "$ROOT/pannexin/project/metadata/species_config.csv" "$OUT/pannexin/project/metadata/" 2>/dev/null || true
if [[ -d "$ROOT/pannexin/project/results/pannexin_path" ]]; then
  "${RSYNC[@]}" "$ROOT/pannexin/project/results/pannexin_path" "$OUT/pannexin/project/results/"
fi
for panx_site in panx_vs_inx pannexin_similarity pannexin_panel pannexin_phylogeny pannexin_curator pannexin_clade_comparison pannexin_annotation_status pannexin_insights pannexin_exon_notes; do
  if [[ -d "$ROOT/pannexin/project/results/$panx_site" ]]; then
    echo "  + pannexin/results/$panx_site"
    "${RSYNC[@]}" --exclude 'work/' "$ROOT/pannexin/project/results/$panx_site" "$OUT/pannexin/project/results/"
  fi
done
# light copies of reference FASTAs (small UniProt panel)
if [[ -d "$ROOT/pannexin/project/data/references/pannexins" ]]; then
  mkdir -p "$OUT/pannexin/project/data/references"
  "${RSYNC[@]}" "$ROOT/pannexin/project/data/references/pannexins" "$OUT/pannexin/project/data/references/"
fi
touch "$OUT/pannexin/project/data/references/pannexins/.gitkeep" 2>/dev/null || true
touch "$OUT/pannexin/project/data/annotations/pannexin/.gitkeep"

# --- share front door ---
cat > "$OUT/README.md" <<'EOF'
# Gap-junction comparative genomics — GitHub share

Lean copy of the thesis websites and scripts (~tens of MB).  
**Not** included: genomes, SynVoy run folders, per-species discovery dumps (~400 GB on the full machine).

## Open the website (local)

HTML must be opened in a browser (GitHub shows source only).

**Start here (clearest overview — now inx → panx → cnx):**

[`project/results/phylogenetic_story/index.html`](project/results/phylogenetic_story/index.html)

**Main guided path:**

[`project/results/gap_junction_path/index.html`](project/results/gap_junction_path/index.html)

**Recovery of unannotated innexins:**

[`project/results/new_species_gallery/index.html`](project/results/new_species_gallery/index.html)

**Innexin × connexin contrast:**

[`project/results/inx_vs_cnx/index.html`](project/results/inx_vs_cnx/index.html)

**Pannexin × innexin:**

[`pannexin/project/results/panx_vs_inx/index.html`](pannexin/project/results/panx_vs_inx/index.html)

**Pannexin path / insights / clade / status / similarity / phylogeny / gaps:**

[`pannexin/project/results/pannexin_path/index.html`](pannexin/project/results/pannexin_path/index.html) ·
[`insights`](pannexin/project/results/pannexin_insights/index.html) ·
[`clade × type`](pannexin/project/results/pannexin_clade_comparison/index.html) ·
[`annotation status`](pannexin/project/results/pannexin_annotation_status/index.html) ·
[`similarity`](pannexin/project/results/pannexin_similarity/index.html) ·
[`phylogeny`](pannexin/project/results/pannexin_phylogeny/index.html) ·
[`gaps`](pannexin/project/results/pannexin_curator/index.html) ·
[`exon notes`](pannexin/project/results/pannexin_exon_notes/index.html)

**Literature cited:**

[`project/results/sources/index.html`](project/results/sources/index.html)

Full reading order: [`project/results/README.md`](project/results/README.md)

## What is in this folder

| Path | Contents |
|------|----------|
| `project/results/` | HTML result sites and figures |
| `project/results/gene_curator_probe/` | Summary tables, **query pack**, per-species `hits.csv` + `present_products.fasta`, accepted locus coordinates, tool versions |
| `project/data/references/innexins/` | Innexin UniProt reference FASTAs (query sources) |
| `project/data/references/connexins/` | Connexin reference FASTAs |
| `project/metadata/` | Species panels and summary tables |
| `pipeline/` | Python site builders |
| `tools/` | Discovery / download / exon CLIs |
| `scripts/` | Shell runners |
| `tests/` | Pytest |
| `docs/` | Pipeline notes |
| `pannexin/` | Pannexin path site + builders + pannexin refs |

## Reproducibility (curator)

- Query pack: [`project/results/gene_curator_probe/batch_queries.fasta`](project/results/gene_curator_probe/batch_queries.fasta)
- Innexin refs: [`project/data/references/innexins/`](project/data/references/innexins/)
- Accepted locus coordinates: [`project/results/gene_curator_probe/accepted_curator_loci.csv`](project/results/gene_curator_probe/accepted_curator_loci.csv)
- Tool versions: [`project/results/gene_curator_probe/TOOL_VERSIONS.md`](project/results/gene_curator_probe/TOOL_VERSIONS.md)
- Per species: `hits.csv` + `present_products.fasta` under `gene_curator_probe/<Species>/`

## What was left out (on purpose)

- Genome FASTAs and SynVoy run folders (~400 GB)
- `project/results/*_discovery/` heavy per-species search dumps
- Curator `locate_*` work directories / BLAST DBs (regenerable from genomes + query pack)

Identities reported on the sites are **MMseqs2** pairwise identities, not BLAST.

## Rebuild this folder from the full project

```bash
./scripts/build_github_share.sh
# or: ./scripts/build_github_share.sh /path/to/out
```

EOF

# Keep the build script inside the share (paths still make sense when run from full repo)
mkdir -p "$OUT/scripts"
cp "$ROOT/scripts/build_github_share.sh" "$OUT/scripts/"

echo
echo "Done."
du -sh "$OUT"
find "$OUT" -type f | wc -l | awk '{print $1 " files"}'
echo
echo "Push only this folder, e.g.:"
echo "  cd $OUT && git init && git add . && git commit -m 'GitHub share: websites and scripts'"
echo "  # then create a new empty GitHub repo and: git remote add origin ... && git push -u origin main"
