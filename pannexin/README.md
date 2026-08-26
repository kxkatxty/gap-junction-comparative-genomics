# Pannexin comparative genomics

Separate workspace from the innexin / connexin gap-junction project.

**Entry:** [project/results/pannexin_path/index.html](project/results/pannexin_path/index.html)

## Why separate?

Pannexins are chordate / vertebrate homologs of innexins — the third family beside invertebrate innexins and vertebrate connexins. They share the innexin fold (four TM helices, Cys motif 2+2=4), but in vertebrates they act as single-membrane channels (often ATP release), not as classical intercellular gap junctions. Docking is blocked by N-linked glycosylation (NGS) on the extracellular loops. Welzel & Schuster 2022 frame this as the remnant of an early chordate innexin bottleneck; connexins then filled the vertebrate gap-junction niche as a separate sequence family.

This folder keeps its own `pipeline/` and `project/results/` so links and builds do not collide with the gap-junction path.

## Layout

| Folder | Role |
|--------|------|
| `docs/` | Notes and methods |
| `pipeline/` | Python builders (`python3 -m pipeline.<name>`) |
| `scripts/` | Shell runners |
| `project/data/` | References, GFFs, genomes |
| `project/metadata/` | Species panel and queues |
| `project/results/` | HTML sites and tables |
| `tests/` | Pytest |

## Planned chapters

1. What pannexins are (and are not)
2. Three vertebrate paralogs: **PANX1 · PANX2 · PANX3** (teleosts may have a fourth)
3. Sequence similarity within pannexins
4. Homology to **innexins** (not connexins) — [`panx_vs_inx`](project/results/panx_vs_inx/index.html)

## Build

```bash
cd pannexin
./scripts/run_panx_vs_inx.sh          # panel + similarity + panx×inx + path
# or stepwise:
python3 -m pipeline.pannexin_panel_builder
python3 -m pipeline.pannexin_similarity_builder
python3 -m pipeline.panx_vs_inx_builder
./scripts/run_build_path.sh
```

Download UniProt references (from `pannexin/` so paths stay local):

```bash
cd pannexin
python3 ../tools/download_by_family.py project/metadata/species_config.csv
# non-mammal fill (TrEMBL):
python3 ../tools/download_by_family.py project/metadata/species_config_unreviewed_fill.csv
```

Config: `project/metadata/species_config.csv`
