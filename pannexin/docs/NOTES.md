# Pannexin project notes

## Biology snapshot

- **Pannexins** = chordate / vertebrate homologs of **innexins** (third family beside innexins and connexins).
- Typical jawed-vertebrate set: **PANX1, PANX2, PANX3**.
- Usual job in vertebrates: plasma-membrane channels (often ATP), **not** docking gap junctions.
- **Sequence markers:** Cys motif 2+2=4 (inx/panx) vs 3+3=6 (cnx); extracellular **N-X-S/T** as non-docking marker. Stoichiometry 8/6/7 = literature labels only.
- **Annotation issues:** GFF δ-GJD fusion, TSA 1-exon models, UniProt copy-number sampling, Inx2≠inx-2, Xenopus/GJA4 microsynteny — see parent gap-junction path `#artifacts`.
- **Docking blockade (literature):** N-linked glycosylation (NGS) on extracellular loops.
- **Bottleneck (literature):** Welzel & Schuster 2022 (*eLife*) — early chordates lost innexin diversity; lancelets have one glycosylated innexin and no connexins; vertebrate pannexins = retained glycosylated lineage (PANX1–3 after WGDs).
- **Connexins** are a separate sequence family that *do* form vertebrate gap junctions.

## Separation rule

Do not write pannexin HTML into `../project/results/`. Keep everything under `pannexin/project/results/`.

Optional later (deferred): SynVoy neighbourhoods / genome probes — do not run from light rebuilds.

## Resource limits (laptop-safe)

- max **1 thread** by default (`PANX_THREADS`)
- phylogeny **reuses** existing treefile (no recompute in normal rebuilds)
- **no genome miniprot** in the curator page
- **SynVoy is deferred** — not wired into pannexin rebuild scripts

```bash
cd pannexin
./scripts/run_panx_expand.sh   # light only
```

## Modules

1. `pannexin_panel_builder` — UniProt FASTA → curated panel
2. `pannexin_similarity_builder` — within-family MMseqs
3. `panx_vs_inx_builder` — vs innexins only
4. `pannexin_phylogeny_builder` — reuse existing tree
5. `pannexin_curator_builder` — annotation-gap status table only
6. `pannexin_path_builder` — guided HTML path

**Later (not now):** SynVoy neighbourhoods, skate genome miniprot, heavier IQ-TREE.
