# Gene curator probe — innexin species with zero accepted pipeline hits

**Skill used:** `gene-curation-skills-main/skills/curating-gene-families`  
**Mode:** Stage 1–2–6 micro-run (single-genome re-probe; Stages 3–5 panel synteny skipped)  
**Date:** 2026-08-09

## What the skill does (short)

`curating-gene-families` is **not a discovery pipeline**. It is a **genome-direct curation workflow**:

1. Treat NCBI/GFF/miniprot/mmseqs labels as **hypotheses**, not truth  
2. Locate candidates with miniprot only as a **coarse pointer**  
3. **faidx + translate** the locus yourself and eyeball the protein product  
4. File a structural verdict: **present / degraded / fragmentary / absent**

Pipelines harvest; the curator **verifies or recovers** models the pipeline discarded.

## Prior-work check

| Species | Pipeline status | Prior artifacts |
|---------|-----------------|-----------------|
| *Rotaria macrura* | 3 candidates, **0 accepted** (all `too_short`) | `innexin_discovery/Rotaria_macrura/` |
| *Brachionus manjavacas* | 1 candidate, **0 accepted** | `innexin_discovery/Brachionus_manjavacas/` |
| *Adineta vaga* (positive control) | **2 accepted** (weak_manual_review) | same folder family |

Expectation: bdelloid / monogonont rotifers should carry innexins (sister *Adineta* already has hits).

## Method this probe used

- **Queries:** Adineta accepted candidates + Dmel Inx2/shakB + Cele inx-2/unc-9  
- **Trace:** `miniprot --outs=0.5` against genome FASTA  
- **Product:** translate miniprot CDS exons genome-direct; score length / Cys / stop-free ORF  
- **Outputs:** `project/results/gene_curator_probe/{Rotaria_macrura,Brachionus_manjavacas}/`

## Results — *Rotaria macrura*

| Locus | Identity | Product | Curator verdict |
|-------|----------|---------|-----------------|
| OENT01000045.1:119094–122499 | 0.36 vs UNC-9 | **386 aa**, 0 stop, 8 Cys | **PRESENT (missed by filter)** |
| OENT01000325.1:81075–83639 | 0.36 vs UNC-9 | **387 aa**, 0 stop, 7 Cys | **PRESENT (missed by filter)** |
| OENT01001445.1:20220–22805 | 0.36 vs UNC-9 | **370 aa**, 0 stop, 8 Cys | **PRESENT (missed by filter)** |
| OENT01000099.1:209084–209674 | 0.88 vs Adineta short | 197 aa | **FRAGMENTARY**; window ORF ~307 aa recoverable |
| OENT01000318.1:41780–42364 | 0.81 vs Adineta short | 195 aa | **FRAGMENTARY** |

**Why the pipeline said zero:** it kept only short models (104–158 aa) and rejected them as `too_short` / `few_tm_helices`.  
**tblastn had already hit** the same contigs (homology was present); the **gene-model / length filter** dropped the real products.

## Results — *Brachionus manjavacas*

Miniprot recovered **7 loci**, including UNC-9-like hits at ~0.32–0.42 identity. Genome-direct products ≥200 aa written to `genome_direct_products.fasta`. Pipeline had only a 118 aa rejected fragment — same pattern as Rotaria.

## Conclusion

1. **Rotaria / Brachionus are not innexin-absent** under curator re-probe — they look **present but under-modeled**.  
2. The failure mode is **filter / incomplete gene models**, not “no homology in the genome”.  
3. Next step: re-annotate the UNC-9-like loci with genome-direct exon walking (GT–AG), then re-rank; do **not** re-run the whole discovery pipeline from scratch without changing length/TM gates for divergent rotifer innexins.

## Files

- `project/results/gene_curator_probe/Rotaria_macrura/miniprot.gff`  
- `project/results/gene_curator_probe/Rotaria_macrura/genome_direct_products.fasta`  
- `project/results/gene_curator_probe/Brachionus_manjavacas/miniprot.gff`  
- `project/results/gene_curator_probe/Brachionus_manjavacas/genome_direct_products.fasta`  
- This report: `project/results/gene_curator_probe/CURATOR_PROBE_REPORT.md`
