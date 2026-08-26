# Shared flanking genes — conclusion (GFF microsynteny)

## Main finding

Vertebrate **connexin** loci show **repetitive, conserved neighbor genes** across species (same gene symbols in multiple genomes). **Innexin** loci do **not**: insects and nematodes mostly use species-specific labels (CG*, LOC*), so the same neighbors rarely appear under identical names across species.

## Connexins

- Panels analysed: **24**
- High conservation (≥3 neighbors in all species): **16**
- Moderate: **4** · Low: **0** · None: **4**

**Strongest blocks** (neighbors shared in every species of the panel):

- **GJB1** (4 species): foxo4, il2rg, itgb1bp2, med12, nlgn3, nono, snx12, taf1, zmym3
- **Gjb1** (2 species): foxo4, il2rg, itgb1bp2, med12, nlgn3, nono, snx12, taf1, zmym3
- **Gjc1** (2 species): adam11, c1ql1, ccdc103, ccdc43, eftud2, fam187a, gfap, higd1b, kif18b
- **GJA4** (3 species): dlgap3, gjb3, gjb4, gjb5, smim12
- **GJB3** (2 species): dlgap3, gja4, gjb4, gjb5, smim12
- **Gja4** (2 species): dlgap3, gjb3, gjb4, gjb5, smim12
- **Gjb3** (2 species): dlgap3, gja4, gjb4, gjb5, smim12
- **Gjb4** (2 species): dlgap3, gja4, gjb3, gjb5, smim12

Recurrent connexin-neighbor patterns:
- **GJA1** locus: **TBC1D32** appears beside GJA1 in all mammal/bird panels.
- **GJA3 / GJB6 / Gja3**: shared **GJB2, GJB6, ZMYM2** block (connexin cluster).
- **GJA4 / GJB3 / GJB4 / GJB5**: shared **DLGAP3, GJB3–GJB5, SMIM12** — classic mammalian connexin β-cluster on one chromosome.
- **GJB1** (X-linked): **FOXO4, IL2RG, NLGN3, MED12** neighbors repeat across primates/rodents.

## Innexins

- Panels analysed: **2**
- High / moderate conservation: **0**

- **inx_19** (2 species): 0 neighbors in all species; best partial overlap = none (tier: none).
- **shakB** (3 species): 0 neighbors in all species; best partial overlap = none (tier: none).

## Interpretation

1. **Connexin microsynteny is conserved** in vertebrates — the same flanking genes (often other connexins or TBC1D32) recur across ortholog panels.
2. **Innexin microsynteny is not symbol-conserved** across insects/nematodes — different annotation systems hide homology; local duplication (Dmel cluster) matters more than cross-species neighbour identity.
3. **Gene-order plots** should be read clade-by-clade: shared neighbour names are informative for connexins; for innexins, compare architecture and SynVoy case studies rather than gene symbols alone.
