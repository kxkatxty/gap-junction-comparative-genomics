# Case Study 2: Ortholog synteny across related genomes

Compares **Inx2** and **shakB** from *Drosophila melanogaster* against SynVoy-selected related genomes.

## Target species searched
- *Drosophila helvetica*, *D. phalerata*, *D. busckii*, *D. limbata* (Drosophilidae)
- *Anopheles rivulorum* (Diptera; more distant; Inx2 only)

## Support classes
- `confident_ortholog` — high-confidence GOI
- `probable_ortholog` — probable GOI
- `weak_goi_signal` — GOI candidates, ambiguous only
- `flanking_synteny_only` — conserved neighborhood, no GOI

## Summary
- **Inx2**: searched 5 genomes; best support in *Drosophila limbata* (probable_ortholog)
- **shakB**: searched 3 genomes; best support in *Drosophila busckii* (confident_ortholog)

## Tables
- `case2_ortholog_support.csv`
- `case2_gene_summary.csv`
- `case2_inx2_vs_shakB.csv`

## Figures
Multi-species synteny HTML plots in `figures/`.
