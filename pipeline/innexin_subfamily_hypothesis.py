#!/usr/bin/env python3
"""Four innexin subfamilies from the reference-panel tree, plus flanking-gene comparison.

Suggested grouping (from innexin_reference_panel tree):
  SF1  shakB clade
  SF2  Inx1 / ogre clade
  SF3  Inx3 + Inx7 + all nematode innexins
  SF4  Inx2 + Inx4(zpg) + Inx5 + Inx6 (insect expansion)

Use these clades for synteny comparisons (not shared names like Inx2 vs inx-2).
"""

from __future__ import annotations

import argparse
import csv
import re
from collections import defaultdict
from pathlib import Path

from pipeline.common import RESULTS_DIR, write_csv
from pipeline.foxp2_style_synteny import (
    LocusSpec,
    _innexin_gff,
    build_panel,
    load_row,
)

OUT = RESULTS_DIR / "innexin_subfamilies"
PHYLO = RESULTS_DIR / "synteny_phylo_guided"
FOX = RESULTS_DIR / "synteny_foxp2_style"

# Tip → subfamily (Dmel + Cele anchors)
SUBFAMILIES: dict[str, dict[str, str]] = {
    "SF1_shakB": {
        "label": "shakB clade",
        "dmel_genes": "shakB",
        "members": "Dmel shakB; mosquito shakB orthologs",
        "genomic_note": "Distal X cluster ~20.76 Mb (NC_004354.4)",
    },
    "SF2_Inx1_ogre": {
        "label": "Inx1 / ogre clade",
        "dmel_genes": "ogre (Inx1)",
        "members": "Dmel ogre; Schistocerca inx1",
        "genomic_note": "Proximal X cluster with Inx7/Inx2 (~6.97 Mb) — genomic neighbor, not SF3 sister",
    },
    "SF3_Inx3_Inx7_nematode": {
        "label": "Inx3 + Inx7 + nematode innexins",
        "dmel_genes": "Inx3, Inx7",
        "members": "Dmel Inx3, Inx7; all C. elegans innexins (inx-*, unc-7/9, eat-5)",
        "genomic_note": "Inx7 sits in proximal X; Inx3 on 3R — tree groups them with nematodes despite split loci",
    },
    "SF4_Inx2_expansion": {
        "label": "Inx2 / Inx4–6 insect expansion",
        "dmel_genes": "Inx2, zpg(Inx4), Inx5, Inx6",
        "members": "Dmel Inx2/4/5/6; Schistocerca inx2; Anopheles ZPG/Inx4",
        "genomic_note": "Inx2 in proximal X tandem with ogre/Inx7; Inx4–6 elsewhere",
    },
}

ASSIGNMENTS = [
    # tip_key, gene, organism, subfamily
    ("INX1_DROME", "ogre", "Drosophila melanogaster", "SF2_Inx1_ogre"),
    ("INX2_DROME", "Inx2", "Drosophila melanogaster", "SF4_Inx2_expansion"),
    ("INX3_DROME", "Inx3", "Drosophila melanogaster", "SF3_Inx3_Inx7_nematode"),
    ("INX4_DROME", "zpg", "Drosophila melanogaster", "SF4_Inx2_expansion"),
    ("INX5_DROME", "Inx5", "Drosophila melanogaster", "SF4_Inx2_expansion"),
    ("INX6_DROME", "Inx6", "Drosophila melanogaster", "SF4_Inx2_expansion"),
    ("INX7_DROME", "Inx7", "Drosophila melanogaster", "SF3_Inx3_Inx7_nematode"),
    ("SHAKB_DROME", "shakB", "Drosophila melanogaster", "SF1_shakB"),
    ("INX2_CAEEL", "inx-2", "Caenorhabditis elegans", "SF3_Inx3_Inx7_nematode"),
    ("INX3_CAEEL", "inx-3", "Caenorhabditis elegans", "SF3_Inx3_Inx7_nematode"),
    ("UNC9_CAEEL", "unc-9", "Caenorhabditis elegans", "SF3_Inx3_Inx7_nematode"),
    ("UNC7_CAEEL", "unc-7", "Caenorhabditis elegans", "SF3_Inx3_Inx7_nematode"),
]


def _family_stem(name: str) -> str:
    """Collapse paralog-ish names to a coarse family stem (ABC1/ABC2 → ABC)."""
    n = name.strip()
    if not n or n.startswith("lncRNA") or n.startswith("mir-") or n.startswith("snoRNA"):
        return ""
    # CG12345 → keep CG stem for FlyBase anonymous genes
    if re.fullmatch(r"CG\d+", n):
        return n
    # strip trailing digits / isoform letters: ABC1, ABC-2, gjb3 → stem
    stem = re.sub(r"[-_]?\d+[A-Za-z]?$", "", n)
    stem = re.sub(r"\.(L|S|a|b|c)$", "", stem)
    return stem or n


def write_hypothesis_tables(out: Path) -> None:
    out.mkdir(parents=True, exist_ok=True)
    write_csv(
        out / "subfamily_definitions.csv",
        [
            {
                "subfamily_id": k,
                "label": v["label"],
                "dmel_genes": v["dmel_genes"],
                "members": v["members"],
                "genomic_note": v["genomic_note"],
            }
            for k, v in SUBFAMILIES.items()
        ],
    )
    write_csv(
        out / "tip_subfamily_assignments.csv",
        [
            {
                "tip_key": a[0],
                "gene": a[1],
                "organism": a[2],
                "subfamily_id": a[3],
                "subfamily_label": SUBFAMILIES[a[3]]["label"],
            }
            for a in ASSIGNMENTS
        ],
    )


def flanking_family_analysis(out: Path) -> Path:
    """Compare flanking gene-name stems around the two Dmel X clusters + Inx3."""
    dmel = _innexin_gff("Drosophila_melanogaster")
    if not dmel:
        raise FileNotFoundError("Dmel GFF missing")

    loci = [
        (
            "proximal_X_SF_mixed",
            "Proximal X (ogre SF2 · Inx7 SF3 · Inx2 SF4)",
            LocusSpec(
                "proximal",
                dmel,
                "NC_004354.4",
                ("gene-Dmel_CG3039", "ogre", "gene-Dmel_CG2977", "Inx7", "gene-Dmel_CG4590", "Inx2"),
                6_840_000,
                7_100_000,
                clade="insect",
            ),
        ),
        (
            "distal_X_SF1",
            "Distal X (shakB SF1)",
            LocusSpec(
                "distal",
                dmel,
                "NC_004354.4",
                ("gene-Dmel_CG34358", "shakB"),
                20_700_000,
                20_950_000,
                clade="insect",
            ),
        ),
        (
            "Inx3_3R_SF3",
            "Inx3 on 3R (SF3)",
            LocusSpec(
                "inx3",
                dmel,
                "NT_033777.3",
                ("gene-Dmel_CG1448", "Inx3"),
                28_820_000,
                28_900_000,
                clade="insect",
            ),
        ),
    ]

    per_locus: dict[str, set[str]] = {}
    detail_rows: list[dict[str, str]] = []
    for lid, label, spec in loci:
        _, genes = load_row(spec)
        stems: set[str] = set()
        for g in genes:
            if g.is_goi:
                detail_rows.append(
                    {
                        "locus_id": lid,
                        "locus_label": label,
                        "gene_name": g.name,
                        "family_stem": g.name,
                        "role": "innexin_goi",
                        "start": str(g.start),
                        "end": str(g.end),
                    }
                )
                continue
            stem = _family_stem(g.name)
            if not stem:
                continue
            stems.add(stem)
            detail_rows.append(
                {
                    "locus_id": lid,
                    "locus_label": label,
                    "gene_name": g.name,
                    "family_stem": stem,
                    "role": "flank",
                    "start": str(g.start),
                    "end": str(g.end),
                }
            )
        per_locus[lid] = stems

    # Pairwise shared family stems between loci (WGD / segmental duplication signal)
    ids = list(per_locus.keys())
    share_rows: list[dict[str, str]] = []
    for i, a in enumerate(ids):
        for b in ids[i + 1 :]:
            shared = sorted(per_locus[a] & per_locus[b])
            share_rows.append(
                {
                    "locus_a": a,
                    "locus_b": b,
                    "n_shared_family_stems": str(len(shared)),
                    "n_a": str(len(per_locus[a])),
                    "n_b": str(len(per_locus[b])),
                    "shared_stems": ";".join(shared[:40]) + ("..." if len(shared) > 40 else ""),
                }
            )

    write_csv(out / "dmel_flanking_genes_detail.csv", detail_rows)
    write_csv(out / "dmel_flanking_family_sharing.csv", share_rows)
    return out / "dmel_flanking_family_sharing.csv"


def build_sf3_vs_nematode_panel(out: Path) -> None:
    """Tree-guided: Inx7 (proximal) + Inx3 (3R) vs nematode SF3 members."""
    dmel = _innexin_gff("Drosophila_melanogaster")
    cele = _innexin_gff("Caenorhabditis_elegans")
    if not dmel or not cele:
        return
    specs = [
        LocusSpec(
            "Dmel Inx7 (SF3) — proximal X",
            dmel,
            "NC_004354.4",
            ("gene-Dmel_CG2977", "Inx7"),
            6_980_000,
            7_010_000,
            clade="insect",
        ),
        LocusSpec(
            "Dmel Inx3 (SF3) — 3R",
            dmel,
            "NT_033777.3",
            ("gene-Dmel_CG1448", "Inx3"),
            28_848_000,
            28_865_000,
            clade="insect",
        ),
        LocusSpec(
            "Cele inx-2 (SF3)",
            cele,
            "NC_003284.9",
            ("gene-CELE_F08G12.10", "inx-2"),
            11_290_000,
            11_320_000,
            clade="nematode",
        ),
        LocusSpec(
            "Cele inx-3 (SF3)",
            cele,
            "NC_003284.9",
            ("gene-CELE_F22F4.2", "inx-3"),
            5_980_000,
            6_020_000,
            clade="nematode",
        ),
        LocusSpec(
            "Cele unc-9 (SF3)",
            cele,
            "NC_003284.9",
            ("gene-CELE_R12H7.1", "unc-9"),
            13_200_000,
            13_230_000,
            clade="nematode",
        ),
        LocusSpec(
            "Cele unc-7 (SF3)",
            cele,
            "NC_003284.9",
            ("gene-CELE_R07D5.1", "unc-7"),
            15_115_000,
            15_155_000,
            clade="nematode",
        ),
    ]
    build_panel(
        "SF3_Inx3_Inx7_vs_nematode_innexins",
        "Tree-guided SF3: Dmel Inx3/Inx7 vs nematode innexins (not name-matched Inx2)",
        specs,
        out,
    )


def write_readme(out: Path, share_path: Path) -> None:
    share = list(csv.DictReader(share_path.open())) if share_path.exists() else []
    lines = [
        "# Innexin subfamilies (from the reference tree)",
        "",
        "Grouping taken from the innexin reference-panel phylogeny",
        "(`project/results/phylogeny/innexin_reference_panel.treefile`).",
        "",
        "## Four subfamilies",
        "",
    ]
    for k, v in SUBFAMILIES.items():
        lines.append(f"### {k} — {v['label']}")
        lines.append(f"- **Dmel genes:** {v['dmel_genes']}")
        lines.append(f"- **Clade members:** {v['members']}")
        lines.append(f"- **Genomic note:** {v['genomic_note']}")
        lines.append("")
    lines.extend(
        [
            "## How this guides synteny",
            "",
            "1. **Do not** treat insect Inx2 and nematode inx-2 as orthologs for synteny.",
            "2. **Do** compare **SF3** (Inx3 + Inx7) neighbourhoods to **nematode innexins** as a group.",
            "3. Drosophila has **two genomic clusters on X** (proximal ogre–Inx7–Inx2; distal shakB)",
            "   plus **Inx3 on 3R** — genomic clustering ≠ the four tree subfamilies one-to-one.",
            "4. Shared flanking **gene families** between clusters would support segmental / WGD-like",
            "   duplication of chromosomal pieces (still to test; see flanking tables).",
            "",
            "## Flanking family sharing (Dmel loci)",
            "",
        ]
    )
    for row in share:
        lines.append(
            f"- `{row['locus_a']}` ↔ `{row['locus_b']}`: "
            f"**{row['n_shared_family_stems']}** shared family stems "
            f"(of {row['n_a']} vs {row['n_b']})"
        )
    lines.extend(
        [
            "",
            "## Files",
            "",
            "- `subfamily_definitions.csv`",
            "- `tip_subfamily_assignments.csv`",
            "- `dmel_flanking_genes_detail.csv`",
            "- `dmel_flanking_family_sharing.csv`",
            "- `SF3_Inx3_Inx7_vs_nematode_innexins.png` — SF3 vs nematode synteny panel",
            "",
            "## Status",
            "",
            "Still provisional — check against domain architecture, expression, and SynVoy",
            "panels as they finish.",
            "",
        ]
    )
    (out / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_index(out: Path) -> None:
    png = out / "SF3_Inx3_Inx7_vs_nematode_innexins.png"
    img = (
        f"<div class='card'><img src='{png.name}' alt='SF3 panel'><p>SF3 vs nematode</p></div>"
        if png.exists()
        else ""
    )
    (out / "index.html").write_text(
        f"""<!DOCTYPE html><html><head><meta charset='utf-8'>
<title>Innexin subfamilies</title>
<style>
body{{font-family:system-ui,sans-serif;margin:2rem;background:#f4f7fa;color:#1a2836}}
.card{{background:#fff;padding:1rem;border-radius:12px;box-shadow:0 2px 12px rgba(0,0,0,.06);margin:1rem 0}}
img{{max-width:100%}} h1{{color:#2E86AB}} a{{color:#1a5276}}
</style></head><body>
<h1>Innexin subfamilies</h1>
<p>Four subfamilies from the reference tree, used for synteny. See <a href='README.md'>README.md</a>.</p>
{img}
<p>Also: <a href='../synteny_phylo_guided/index.html'>phylo-guided gallery</a> ·
<a href='../phylogenetic_story/index.html'>phylogeny story</a></p>
</body></html>""",
        encoding="utf-8",
    )


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", type=Path, default=OUT)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    write_hypothesis_tables(args.outdir)
    share = flanking_family_analysis(args.outdir)
    build_sf3_vs_nematode_panel(args.outdir)
    write_readme(args.outdir, share)
    write_index(args.outdir)
    # Mirror SF3 panel into phylo-guided gallery folder
    PHYLO.mkdir(parents=True, exist_ok=True)
    for ext in (".png", ".html"):
        src = args.outdir / f"SF3_Inx3_Inx7_vs_nematode_innexins{ext}"
        if src.exists():
            (PHYLO / src.name).write_bytes(src.read_bytes())
    print(f"Wrote hypothesis + flanking analysis → {args.outdir}")
    print(f"SF3 panel → {args.outdir / 'SF3_Inx3_Inx7_vs_nematode_innexins.png'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
