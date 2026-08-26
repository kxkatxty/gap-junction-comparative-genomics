# Key findings

- MMseqs2 sweep: highest mean subfamily purity 0.89 at min-seq-id=0.40 (16 clusters).
- Within-subfamily median identity ~46.2% vs between-subfamily median 27.9% — groups are sequence-separable but still one family.
- 3-mer embedding cosine: same-subfamily mean 0.279 vs different-subfamily 0.049.
- AA-composition cosine is weaker at separating groups (same 0.959 vs diff 0.955) — composition alone is not enough; k-mers / MMseqs carry more signal.
- mmseqs_pident: SF3↔SF4 mean similarity 28.4782 (n=1702) — tests insect expansion vs nematode-like clade.
- Note: curator subfamily labels come from best-hit mapping; embedding/MMseqs separation supports the tree groups, not final orthology.
