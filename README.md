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

