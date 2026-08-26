# Innexin subfamilies (from the reference tree)

Grouping taken from the innexin reference-panel phylogeny
(`project/results/phylogeny/innexin_reference_panel.treefile`).

## Four subfamilies

### SF1_shakB — shakB clade
- **Dmel genes:** shakB
- **Clade members:** Dmel shakB; mosquito shakB orthologs
- **Genomic note:** Distal X cluster ~20.76 Mb (NC_004354.4)

### SF2_Inx1_ogre — Inx1 / ogre clade
- **Dmel genes:** ogre (Inx1)
- **Clade members:** Dmel ogre; Schistocerca inx1
- **Genomic note:** Proximal X cluster with Inx7/Inx2 (~6.97 Mb) — genomic neighbor, not SF3 sister

### SF3_Inx3_Inx7_nematode — Inx3 + Inx7 + nematode innexins
- **Dmel genes:** Inx3, Inx7
- **Clade members:** Dmel Inx3, Inx7; all C. elegans innexins (inx-*, unc-7/9, eat-5)
- **Genomic note:** Inx7 sits in proximal X; Inx3 on 3R — tree groups them with nematodes despite split loci

### SF4_Inx2_expansion — Inx2 / Inx4–6 insect expansion
- **Dmel genes:** Inx2, zpg(Inx4), Inx5, Inx6
- **Clade members:** Dmel Inx2/4/5/6; Schistocerca inx2; Anopheles ZPG/Inx4
- **Genomic note:** Inx2 in proximal X tandem with ogre/Inx7; Inx4–6 elsewhere

## How this guides synteny

1. **Do not** treat insect Inx2 and nematode inx-2 as orthologs for synteny.
2. **Do** compare **SF3** (Inx3 + Inx7) neighbourhoods to **nematode innexins** as a group.
3. Drosophila has **two genomic clusters on X** (proximal ogre–Inx7–Inx2; distal shakB)
   plus **Inx3 on 3R** — genomic clustering ≠ the four tree subfamilies one-to-one.
4. Shared flanking **gene families** between clusters would support segmental / WGD-like
   duplication of chromosomal pieces (still to test; see flanking tables).

## Flanking family sharing (Dmel loci)

- `proximal_X_SF_mixed` ↔ `distal_X_SF1`: **0** shared family stems (of 21 vs 8)
- `proximal_X_SF_mixed` ↔ `Inx3_3R_SF3`: **0** shared family stems (of 21 vs 7)
- `distal_X_SF1` ↔ `Inx3_3R_SF3`: **0** shared family stems (of 8 vs 7)

## Files

- `subfamily_definitions.csv`
- `tip_subfamily_assignments.csv`
- `dmel_flanking_genes_detail.csv`
- `dmel_flanking_family_sharing.csv`
- `SF3_Inx3_Inx7_vs_nematode_innexins.png` — SF3 vs nematode synteny panel

## Status

Still provisional — check against domain architecture, expression, and SynVoy
panels as they finish.

