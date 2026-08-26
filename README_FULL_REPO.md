# Gap-junction comparative genomics

**For GitHub:** use the lean folder [`github_share/`](github_share/) (~46 MB: websites + scripts + metadata). Rebuild with `./scripts/build_github_share.sh`. Do not push the full local tree (~400 GB genomes and discovery dumps).

## Result pages

| Page | Path |
|------|------|
| Main path | [`project/results/gap_junction_path/index.html`](project/results/gap_junction_path/index.html) |
| Innexin × connexin | [`project/results/inx_vs_cnx/index.html`](project/results/inx_vs_cnx/index.html) |
| Innexin insights | [`project/results/innexin_insights/index.html`](project/results/innexin_insights/index.html) |
| Connexin insights | [`project/results/connexin_insights/index.html`](project/results/connexin_insights/index.html) |
| Phylogeny story (inx → panx → cnx) | [`project/results/phylogenetic_story/index.html`](project/results/phylogenetic_story/index.html) |
| Pannexin insights | [`pannexin/project/results/pannexin_insights/index.html`](pannexin/project/results/pannexin_insights/index.html) |
| Pannexin clade × type | [`pannexin/project/results/pannexin_clade_comparison/index.html`](pannexin/project/results/pannexin_clade_comparison/index.html) |
| Pannexin annotation status | [`pannexin/project/results/pannexin_annotation_status/index.html`](pannexin/project/results/pannexin_annotation_status/index.html) |
| Pannexin path | [`pannexin/project/results/pannexin_path/index.html`](pannexin/project/results/pannexin_path/index.html) |
| Pannexin × innexin | [`pannexin/project/results/panx_vs_inx/index.html`](pannexin/project/results/panx_vs_inx/index.html) |
| Pannexin similarity | [`pannexin/project/results/pannexin_similarity/index.html`](pannexin/project/results/pannexin_similarity/index.html) |
| Pannexin phylogeny | [`pannexin/project/results/pannexin_phylogeny/index.html`](pannexin/project/results/pannexin_phylogeny/index.html) |
| Pannexin gaps | [`pannexin/project/results/pannexin_curator/index.html`](pannexin/project/results/pannexin_curator/index.html) |
| Sources | [`project/results/sources/index.html`](project/results/sources/index.html) |

Reading order: [`project/results/README.md`](project/results/README.md)

Reported identities are **MMseqs2** pairwise identities, not BLAST.  
Open HTML files in a local browser (GitHub shows source only).

## Repository layout

| Folder | Contents |
|--------|----------|
| `project/results/` | HTML sites, figures, tables |
| `pipeline/` | Python builders |
| `tools/` | Discovery / download / exon CLIs |
| `scripts/` | Shell runners |
| `docs/` | Pipeline notes |
| `project/metadata/` | Species lists and summary tables |
| `pannexin/` | Pannexin workspace |
| `SynVoy/` | Synteny tool (run outputs gitignored) |
| `tests/` | Pytest |

Large local files (genomes, SynVoy run folders) are listed in `.gitignore`.

## GitHub Pages (optional)

Deploy from branch root. Main path URL:

`https://<user>.github.io/<repo>/project/results/gap_junction_path/`

## Commands

```bash
python3 -m pipeline.gap_junction_path_builder
python3 -m pipeline.inx_vs_cnx_builder
python3 tools/discover_innexins.py --family innexin --help

cd pannexin && ./scripts/run_build_path.sh
cd pannexin && ./scripts/run_panx_vs_inx.sh
```
Config: `project/metadata/species_config.csv`  
Pannexin config: `pannexin/project/metadata/species_config.csv`
