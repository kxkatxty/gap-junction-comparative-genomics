#!/usr/bin/env python3
"""Phylo-guided synteny: Dmel innexin clusters vs nematode innexins.

Builds gene-order panels that place the Drosophila proximal
cluster (ogre–Inx7–Inx2) and distal shakB locus next to C. elegans
inx-2 / inx-3 / unc-9 neighbourhoods — using the phylogenetic split
(insect vs nematode clades), not name identity alone.
"""

from __future__ import annotations

import argparse
import shutil
from pathlib import Path

from pipeline.common import RESULTS_DIR, write_csv
from pipeline.foxp2_style_synteny import LocusSpec, build_panel, _innexin_gff

OUT_DIR = RESULTS_DIR / "synteny_phylo_guided"
SHOWCASE_STORY = RESULTS_DIR / "phylogenetic_story" / "figures"
CASE1_DIR = RESULTS_DIR / "synteny_case_studies" / "case1_dmel_innexin_cluster"


def phylo_panels() -> list[tuple[str, str, list[LocusSpec]]]:
    dmel = _innexin_gff("Drosophila_melanogaster")
    cele = _innexin_gff("Caenorhabditis_elegans")
    panels: list[tuple[str, str, list[LocusSpec]]] = []
    if not dmel or not cele:
        return panels

    # Panel A: two Dmel X clusters with gene names (one row each)
    panels.append(
        (
            "dmel_two_clusters_named",
            "Two D. melanogaster innexin clusters on X — gene names labeled",
            [
                LocusSpec(
                    "Proximal X — ogre · Inx7 · Inx2",
                    dmel,
                    "NC_004354.4",
                    ("gene-Dmel_CG3039", "ogre", "gene-Dmel_CG2977", "Inx7", "gene-Dmel_CG4590", "Inx2"),
                    6_960_000,
                    7_020_000,
                    clade="insect",
                ),
                LocusSpec(
                    "Distal X — shakB",
                    dmel,
                    "NC_004354.4",
                    ("gene-Dmel_CG34358", "shakB"),
                    20_740_000,
                    20_800_000,
                    clade="insect",
                ),
            ],
        )
    )

    # Panel B: phylo-guided — proximal insect cluster vs nematode innexins
    panels.append(
        (
            "inx_cluster_vs_nematode",
            "Phylo-guided synteny: Dmel proximal cluster (ogre–Inx7–Inx2) vs nematode innexins",
            [
                LocusSpec(
                    "D. melanogaster — proximal (ogre·Inx7·Inx2)",
                    dmel,
                    "NC_004354.4",
                    ("gene-Dmel_CG3039", "ogre", "gene-Dmel_CG2977", "Inx7", "gene-Dmel_CG4590", "Inx2"),
                    6_960_000,
                    7_020_000,
                    clade="insect",
                ),
                LocusSpec(
                    "C. elegans — inx-2",
                    cele,
                    "NC_003284.9",
                    ("gene-CELE_F08G12.10", "inx-2"),
                    11_290_000,
                    11_320_000,
                    clade="nematode",
                ),
                LocusSpec(
                    "C. elegans — inx-3",
                    cele,
                    "NC_003284.9",
                    ("gene-CELE_F22F4.2", "inx-3"),
                    5_980_000,
                    6_020_000,
                    clade="nematode",
                ),
                LocusSpec(
                    "C. elegans — unc-9",
                    cele,
                    "NC_003284.9",
                    ("gene-CELE_R12H7.1", "unc-9"),
                    13_200_000,
                    13_230_000,
                    clade="nematode",
                ),
            ],
        )
    )

    # Panel C: distal / Inx3 side — shakB + Inx3 (3R) vs nematode
    panels.append(
        (
            "inx3_shakB_vs_nematode",
            "Phylo-guided: Dmel shakB (X) + Inx3 (3R) vs C. elegans inx-3 / unc-7",
            [
                LocusSpec(
                    "D. melanogaster — shakB (X)",
                    dmel,
                    "NC_004354.4",
                    ("gene-Dmel_CG34358", "shakB"),
                    20_740_000,
                    20_800_000,
                    clade="insect",
                ),
                LocusSpec(
                    "D. melanogaster — Inx3 (3R)",
                    dmel,
                    "NT_033777.3",
                    ("gene-Dmel_CG1448", "Inx3"),
                    28_848_000,
                    28_865_000,
                    clade="insect",
                ),
                LocusSpec(
                    "C. elegans — inx-3",
                    cele,
                    "NC_003284.9",
                    ("gene-CELE_F22F4.2", "inx-3"),
                    5_980_000,
                    6_020_000,
                    clade="nematode",
                ),
                LocusSpec(
                    "C. elegans — unc-7",
                    cele,
                    "NC_003284.9",
                    ("gene-CELE_R07D5.1", "unc-7"),
                    15_115_000,
                    15_155_000,
                    clade="nematode",
                ),
            ],
        )
    )
    return panels


def write_two_cluster_summary() -> Path:
    """Compact table used by the revised showcase synteny figure."""
    rows = [
        {
            "cluster_id": "proximal_X",
            "cluster_label": "Proximal X ~6.97–7.00 Mb",
            "chrom": "NC_004354.4",
            "innexin_genes": "ogre · Inx7 · Inx2",
            "shared_flanking": "25",
            "unique_flanking": "0",
            "note": "Tandem cluster; SynVoy Case 1 proximal blocks share all flanks",
        },
        {
            "cluster_id": "distal_X",
            "cluster_label": "Distal X ~20.76–20.93 Mb",
            "chrom": "NC_004354.4",
            "innexin_genes": "shakB",
            "shared_flanking": "27",
            "unique_flanking": "0",
            "note": "shakB locus; Inx3 is annotated on 3R (NT_033777.3), not this X cluster",
        },
    ]
    CASE1_DIR.mkdir(parents=True, exist_ok=True)
    path = CASE1_DIR / "case1_two_clusters.csv"
    write_csv(path, rows)
    return path


def write_readme(out_dir: Path, built: list[dict[str, str]]) -> None:
    lines = [
        "# Phylo-guided innexin synteny",
        "",
        "Compares **Drosophila** innexin genomic clusters to **C. elegans** innexin",
        "neighborhoods using the phylogenetic split (insect vs nematode), not gene-name orthology.",
        "",
        "## Biological framing",
        "",
        "- **Proximal X cluster:** ogre (Inx1) – Inx7 – Inx2 (~6.97–7.00 Mb)",
        "- **Distal X cluster:** shakB (~20.76 Mb)",
        "- **Inx3:** annotated on chromosome **3R** (`NT_033777.3`), separate from the X clusters",
        "- **Nematode anchors:** inx-2, inx-3, unc-9, unc-7 (chromosome X / NC_003284.9)",
        "",
        "## Panels",
        "",
    ]
    for row in built:
        lines.append(f"- `{Path(row['png']).name}` — {row['panel']}")
    lines.extend(
        [
            "",
            "## Note",
            "Microsynteny does **not** conserve across the insect–nematode split;",
            "these panels document that phylogenetic depth breaks gene-order conservation.",
            "",
        ]
    )
    (out_dir / "README.md").write_text("\n".join(lines) + "\n", encoding="utf-8")


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--outdir", type=Path, default=OUT_DIR)
    return p.parse_args()


def main() -> int:
    args = parse_args()
    args.outdir.mkdir(parents=True, exist_ok=True)
    write_two_cluster_summary()
    built: list[dict[str, str]] = []
    for panel_id, title, specs in phylo_panels():
        meta = build_panel(panel_id, title, specs, args.outdir)
        built.append(meta)
        print(f"OK  {panel_id} -> {meta['png']}")

    write_readme(args.outdir, built)

    # Copy key panels into story + case1 figures for easy linking
    CASE1_DIR.mkdir(parents=True, exist_ok=True)
    fig_dir = CASE1_DIR / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)
    for name in ("dmel_two_clusters_named.png", "inx_cluster_vs_nematode.png"):
        src = args.outdir / name
        if src.exists():
            shutil.copy2(src, fig_dir / name)
            if SHOWCASE_STORY.exists():
                # keep story's 03 as cluster figure if we regenerate showcase separately
                pass

    index = args.outdir / "index.html"
    cards = "\n".join(
        f"<div class='card'><a href='{Path(r['html']).name}'>"
        f"<img src='{Path(r['png']).name}' alt='{r['panel']}'>"
        f"<p>{r['panel'].replace('_', ' ')}</p></a></div>"
        for r in built
    )
    index.write_text(
        f"""<!DOCTYPE html><html><head><meta charset='utf-8'>
<title>Phylo-guided innexin synteny</title>
<style>
body{{font-family:system-ui,sans-serif;margin:2rem;background:#f4f7fa;color:#1a2836}}
.grid{{display:grid;gap:1.25rem;grid-template-columns:repeat(auto-fit,minmax(320px,1fr))}}
.card{{background:#fff;border-radius:12px;padding:1rem;box-shadow:0 2px 12px rgba(0,0,0,.06)}}
img{{width:100%;border-radius:8px}}
h1{{color:#2E86AB}}
a{{color:#1a5276;text-decoration:none}}
</style></head><body>
<h1>Phylo-guided innexin synteny</h1>
<p>Dmel clusters vs nematode innexins — see <a href='README.md'>README.md</a></p>
<div class='grid'>{cards}</div>
</body></html>""",
        encoding="utf-8",
    )
    print(f"Gallery: {index}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
