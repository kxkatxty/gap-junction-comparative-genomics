#!/usr/bin/env python3
"""Export small curator reproducibility artifacts (products + coordinate table).

Writes, under gene_curator_probe/:
  - <species>/present_products.fasta   (present + fragmentary translations)
  - accepted_curator_loci.csv          (all present loci coordinates)
  - TOOL_VERSIONS.md                   (if --versions also requested here we only do data)

Does not remapping genomes — re-translates overlapping miniprot GFF models.
"""

from __future__ import annotations

import csv
from pathlib import Path

from pipeline.batch_innexin_curator_probe import (
    GENOME_ROOT,
    OUT_ROOT,
    find_genome_fasta,
    parse_miniprot_models,
    translate_model,
)
from pipeline.common import PROJECT_ROOT, write_csv
from pipeline.rescore_curator_el_cys import _pick_gff


def _overlap(a0: int, a1: int, b0: int, b1: int) -> int:
    return max(0, min(a1, b1) - max(a0, b0))


def export_species_products(species: str) -> tuple[int, list[dict[str, str]]]:
    sp_dir = OUT_ROOT / species
    hits_path = sp_dir / "hits.csv"
    if not hits_path.exists():
        return 0, []
    hits = list(csv.DictReader(hits_path.open(encoding="utf-8", newline="")))
    keep = [h for h in hits if h.get("verdict") in {"present", "fragmentary"}]
    genome = find_genome_fasta(GENOME_ROOT / species)
    gff = _pick_gff(sp_dir)
    models = parse_miniprot_models(gff) if gff and genome else []

    fasta_lines: list[str] = []
    coord_rows: list[dict[str, str]] = []
    n_ok = 0
    for i, h in enumerate(keep, 1):
        seqid, start, end = h["seqid"], int(h["start"]), int(h["end"])
        prot = ""
        if genome and models:
            best_m, best_ov = None, 0
            for m in models:
                if m["seqid"] != seqid:
                    continue
                ov = _overlap(start, end, int(m["start"]), int(m["end"]))
                if ov > best_ov:
                    best_ov, best_m = ov, m
            if best_m is not None and best_ov >= 50:
                try:
                    prot = translate_model(genome, best_m)
                except Exception:
                    prot = ""
        clean = (prot[:-1] if prot.endswith("*") else prot).split("*")[0] if prot else ""
        header = (
            f">{species}|{h.get('verdict')}|{seqid}:{start}-{end}|"
            f"{h.get('strand', '+')}|id={h.get('identity','')}|aa={h.get('aa_length','')}|"
            f"el_motif={h.get('el_motif','')}|target={((h.get('target') or '').split() or [''])[0]}"
        )
        if clean:
            fasta_lines.append(f">{header}\n{clean}\n")
            n_ok += 1
        if h.get("verdict") == "present":
            coord_rows.append(
                {
                    "species": species,
                    "seqid": seqid,
                    "start": str(start),
                    "end": str(end),
                    "strand": h.get("strand") or "+",
                    "identity": h.get("identity") or "",
                    "aa_length": h.get("aa_length") or "",
                    "cys": h.get("cys") or "",
                    "el1_cys": h.get("el1_cys") or "",
                    "el2_cys": h.get("el2_cys") or "",
                    "el_motif": h.get("el_motif") or "",
                    "tm_pred": h.get("tm_pred") or "",
                    "target": h.get("target") or "",
                    "product_exported": "yes" if clean else "no",
                }
            )
    out_fa = sp_dir / "present_products.fasta"
    if fasta_lines:
        out_fa.write_text("".join(fasta_lines), encoding="utf-8")
    elif out_fa.exists():
        out_fa.unlink()
    return n_ok, coord_rows


def main() -> int:
    species_dirs = sorted(
        p.name for p in OUT_ROOT.iterdir() if p.is_dir() and (p / "hits.csv").exists()
    )
    all_coords: list[dict[str, str]] = []
    total_fa = 0
    for sp in species_dirs:
        n, rows = export_species_products(sp)
        total_fa += n
        all_coords.extend(rows)
        if rows:
            print(f"{sp}: present={len(rows)} products_written={n}")
    write_csv(OUT_ROOT / "accepted_curator_loci.csv", all_coords)
    print(f"Wrote {OUT_ROOT / 'accepted_curator_loci.csv'} ({len(all_coords)} present loci)")
    print(f"Product sequences exported: {total_fa}")
    print(f"Relative to project: {(OUT_ROOT / 'accepted_curator_loci.csv').relative_to(PROJECT_ROOT)}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
