#!/usr/bin/env python3
"""Drosophila melanogaster innexin recovery benchmark (ground truth = 8 known genes).

Compares:
  A) current curator locate: whole-genome miniprot (floor)
  B) improved locate: MMseqs2 region list → per-region miniprot + guards

Writes recovery tables and a short markdown report for the thesis.
"""

from __future__ import annotations

import argparse
import csv
import json
import re
from dataclasses import dataclass
from pathlib import Path

from pipeline.batch_innexin_curator_probe import build_query_fasta
from pipeline.common import PROJECT_ROOT, RESULTS_DIR
from pipeline.innexin_locate import (
    LocatedLocus,
    default_curator_queries,
    locate_miniprot_only,
    locate_with_mmseqs_then_miniprot,
    present_loci,
)

OUT = RESULTS_DIR / "dmel_innexin_benchmark"
GENOME = (
    PROJECT_ROOT
    / "project/data/genomes/innexin_discovery/Drosophila_melanogaster"
    / "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.fna.decompressed.fna"
)
GFF = (
    PROJECT_ROOT
    / "project/data/annotations/innexin/Drosophila_melanogaster"
    / "GCF_000001215.4_Release_6_plus_ISO1_MT_genomic.gff"
)

# UniProt accessions used in reference FASTA filenames
GENE_TO_ACCESSION = {
    "ogre": "P27716",
    "Inx2": "Q9V427",
    "Inx3": "Q9VAS7",
    "zpg": "Q9VRX6",
    "Inx5": "Q9VWL5",
    "Inx6": "Q9VR82",
    "Inx7": "Q9V3W6",
    "shakB": "P33085",
}


@dataclass
class TruthGene:
    name: str
    seqid: str
    start: int
    end: int
    strand: str
    accession: str


def load_truth_genes(gff: Path = GFF) -> list[TruthGene]:
    aliases = {
        "ogre": {"ogre"},
        "Inx2": {"Inx2"},
        "Inx3": {"Inx3"},
        "zpg": {"zpg"},
        "Inx5": {"Inx5"},
        "Inx6": {"Inx6"},
        "Inx7": {"Inx7"},
        "shakB": {"shakB"},
    }
    found: dict[str, TruthGene] = {}
    with gff.open(encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            m = re.search(r"Name=([^;]+)", parts[8])
            if not m:
                continue
            gname = m.group(1)
            for gene, names in aliases.items():
                if gname in names and gene not in found:
                    found[gene] = TruthGene(
                        name=gene,
                        seqid=parts[0],
                        start=int(parts[3]),
                        end=int(parts[4]),
                        strand=parts[6],
                        accession=GENE_TO_ACCESSION[gene],
                    )
    missing = [g for g in GENE_TO_ACCESSION if g not in found]
    if missing:
        raise RuntimeError(f"Missing truth genes in GFF: {missing}")
    return [found[g] for g in GENE_TO_ACCESSION]


def intervals_overlap(a0: int, a1: int, b0: int, b1: int, *, pad: int = 0) -> bool:
    return not (a1 + pad < b0 or b1 + pad < a0)


def score_recovery(
    truth: list[TruthGene],
    loci: list[LocatedLocus],
    *,
    pad_bp: int = 2000,
) -> list[dict[str, str]]:
    candidates = [x for x in loci if x.verdict in {"present", "fragmentary"}]
    rows: list[dict[str, str]] = []
    for gene in truth:
        hits = [
            loc
            for loc in candidates
            if loc.seqid == gene.seqid
            and intervals_overlap(loc.start, loc.end, gene.start, gene.end, pad=pad_bp)
        ]
        hits.sort(key=lambda x: (-(x.verdict == "present"), -x.identity, -x.aa_length))
        best = hits[0] if hits else None
        rows.append(
            {
                "gene": gene.name,
                "accession": gene.accession,
                "truth_seqid": gene.seqid,
                "truth_start": str(gene.start),
                "truth_end": str(gene.end),
                "recovered": "yes" if best and best.verdict == "present" else "no",
                "best_verdict": best.verdict if best else "",
                "best_identity": f"{best.identity:.4f}" if best else "",
                "best_aa": str(best.aa_length) if best else "",
                "best_exons": str(best.exon_count) if best else "",
                "best_locus": (
                    f"{best.seqid}:{best.start}-{best.end}" if best else ""
                ),
                "n_overlapping_loci": str(len(hits)),
            }
        )
    return rows


def write_csv(path: Path, rows: list[dict[str, str]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        path.write_text("", encoding="utf-8")
        return
    with path.open("w", encoding="utf-8", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()))
        w.writeheader()
        w.writerows(rows)


def write_report(
    path: Path,
    *,
    truth: list[TruthGene],
    old_rows: list[dict[str, str]],
    new_rows: list[dict[str, str]],
    old_present: int,
    new_present: int,
    old_n_loci: int,
    new_n_loci: int,
) -> None:
    old_hit = sum(1 for r in old_rows if r["recovered"] == "yes")
    new_hit = sum(1 for r in new_rows if r["recovered"] == "yes")
    missed_old = [r["gene"] for r in old_rows if r["recovered"] != "yes"]
    missed_new = [r["gene"] for r in new_rows if r["recovered"] != "yes"]
    gained = [
        r["gene"]
        for r in new_rows
        if r["recovered"] == "yes"
        and next(x for x in old_rows if x["gene"] == r["gene"])["recovered"] != "yes"
    ]
    lines = [
        "# Drosophila innexin recovery benchmark",
        "",
        "Ground truth: **8** annotated *D. melanogaster* innexins "
        "(ogre/Inx1, Inx2, Inx3, zpg/Inx4, Inx5, Inx6, Inx7, shakB) "
        "on Release 6 plus ISO1 MT (`GCF_000001215.4`).",
        "",
        "Queries for both methods: the **standard curator query pack** "
        "(not the 8 Dmel self-proteins), matching how non-model genomes are searched.",
        "",
        "## Calibration (thesis number)",
        "",
        f"| Method | Known genes recovered | Present loci called |",
        f"|--------|----------------------:|--------------------:|",
        f"| A. Whole-genome miniprot (current) | **{old_hit}/8** | {old_present} "
        f"(raw kept loci {old_n_loci}) |",
        f"| B. MMseqs2 regions → per-region miniprot | **{new_hit}/8** | {new_present} "
        f"(raw kept loci {new_n_loci}) |",
        "",
        f"Missed by A: {', '.join(missed_old) if missed_old else 'none'}",
        "",
        f"Missed by B: {', '.join(missed_new) if missed_new else 'none'}",
        "",
        f"Gained by B over A: {', '.join(gained) if gained else 'none'}",
        "",
        "## Interpretation",
        "",
        "Method A returns roughly one best model per query protein, so per-species "
        "`present_loci` counts are a **floor**, not a census. On Drosophila, where "
        f"the true answer is 8, method A recovers **{old_hit}/8** "
        f"(short by a factor of roughly {8 / max(old_hit, 1):.1f}).",
        "",
        "Method B first lists homologous genomic windows with MMseqs2 "
        "(translated search), then runs miniprot inside each window with "
        "exon-count / length / chimera guards, so paralogs at distinct loci "
        "are no longer collapsed by the single best-per-query behaviour.",
        "",
        "## Truth set",
        "",
    ]
    for g in truth:
        lines.append(
            f"- **{g.name}** (`{g.accession}`): `{g.seqid}:{g.start}-{g.end}` ({g.strand})"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `truth_genes.csv`",
            "- `method_A_miniprot_only/loci.csv`",
            "- `method_A_recovery.csv`",
            "- `method_B_mmseqs_regions/loci.csv`",
            "- `method_B_mmseqs_regions/regions.csv`",
            "- `method_B_recovery.csv`",
            "- `benchmark_summary.json`",
            "",
        ]
    )
    path.write_text("\n".join(lines), encoding="utf-8")


def main() -> int:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--outs", type=float, default=0.5)
    p.add_argument("--evalue", type=float, default=1e-5)
    p.add_argument(
        "--skip-b",
        action="store_true",
        help="Only run method A (current floor)",
    )
    args = p.parse_args()

    if not GENOME.exists():
        raise SystemExit(f"Missing Dmel genome: {GENOME}")
    if not GFF.exists():
        raise SystemExit(f"Missing Dmel GFF: {GFF}")

    OUT.mkdir(parents=True, exist_ok=True)
    truth = load_truth_genes()
    write_csv(
        OUT / "truth_genes.csv",
        [
            {
                "gene": g.name,
                "accession": g.accession,
                "seqid": g.seqid,
                "start": str(g.start),
                "end": str(g.end),
                "strand": g.strand,
            }
            for g in truth
        ],
    )

    queries = default_curator_queries()
    if not queries.exists() or queries.stat().st_size == 0:
        queries = build_query_fasta(queries)

    print(f"Queries: {queries}")
    print(f"Genome:  {GENOME}")
    print("Running method A: whole-genome miniprot ...")
    old_loci = locate_miniprot_only(
        GENOME,
        queries,
        OUT / "method_A_miniprot_only",
        threads=args.threads,
        outs=args.outs,
    )
    old_rows = score_recovery(truth, old_loci)
    write_csv(OUT / "method_A_recovery.csv", old_rows)
    old_present = len(present_loci(old_loci))
    old_hit = sum(1 for r in old_rows if r["recovered"] == "yes")
    print(f"  Method A: recovered {old_hit}/8 known; present_loci={old_present}")

    new_rows: list[dict[str, str]] = []
    new_loci: list[LocatedLocus] = []
    new_present = 0
    new_hit = 0
    if not args.skip_b:
        print("Running method B: MMseqs2 regions → per-region miniprot ...")
        new_loci = locate_with_mmseqs_then_miniprot(
            GENOME,
            queries,
            OUT / "method_B_mmseqs_regions",
            threads=args.threads,
            evalue=args.evalue,
            outs=args.outs,
        )
        new_rows = score_recovery(truth, new_loci)
        write_csv(OUT / "method_B_recovery.csv", new_rows)
        new_present = len(present_loci(new_loci))
        new_hit = sum(1 for r in new_rows if r["recovered"] == "yes")
        print(f"  Method B: recovered {new_hit}/8 known; present_loci={new_present}")

    summary = {
        "genome": str(GENOME),
        "queries": str(queries),
        "truth_n": len(truth),
        "method_A": {
            "recovered": old_hit,
            "present_loci": old_present,
            "kept_loci": len(old_loci),
            "genes": {r["gene"]: r["recovered"] for r in old_rows},
        },
        "method_B": {
            "recovered": new_hit,
            "present_loci": new_present,
            "kept_loci": len(new_loci),
            "genes": {r["gene"]: r["recovered"] for r in new_rows} if new_rows else {},
        },
    }
    (OUT / "benchmark_summary.json").write_text(
        json.dumps(summary, indent=2) + "\n", encoding="utf-8"
    )
    write_report(
        OUT / "BENCHMARK_REPORT.md",
        truth=truth,
        old_rows=old_rows,
        new_rows=new_rows or old_rows,
        old_present=old_present,
        new_present=new_present,
        old_n_loci=len(old_loci),
        new_n_loci=len(new_loci),
    )
    print(f"Wrote {OUT / 'BENCHMARK_REPORT.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
