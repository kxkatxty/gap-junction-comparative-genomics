#!/usr/bin/env python3
"""Microsynteny gene-order plots from NCBI GFF annotations.

Layout inspired by synteny_foxp2__1_.pdf: each species is one row, genes in genomic
order left-to-right, gene of interest highlighted. Uses annotated GFF windows only
(no SynVoy). Method follows standard microsynteny practice (gene order in locus
blocks; cf. MCScanX / Synteny Portal gene-order views).
"""

from __future__ import annotations

import argparse
import re
import sys
from dataclasses import dataclass
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch

from pipeline.common import ANNOTATIONS_DIR, PROJECT_ROOT, RESULTS_DIR, write_csv

OUT_DIR = RESULTS_DIR / "synteny_foxp2_style"

GOI_COLOR = "#C0392B"
FLANK_COLOR = "#5DADE2"
FLANK_ALT = "#AED6F1"
ROW_LABEL_W = 0.22

# Distinct palette for shared flanking genes (same name → same color across rows).
GENE_PALETTE = (
    "#5DADE2",
    "#58D68D",
    "#F4D03F",
    "#AF7AC5",
    "#48C9B0",
    "#F5B041",
    "#85C1E9",
    "#82E0AA",
    "#F1948A",
    "#D7BDE2",
    "#76D7C4",
    "#F8C471",
    "#7FB3D5",
    "#A9DFBF",
    "#D2B4DE",
    "#F9E79F",
)


@dataclass(frozen=True)
class LocusSpec:
    species_label: str
    gff_path: Path
    chrom: str
    goi_gene_ids: tuple[str, ...]
    window_start: int | None = None
    window_end: int | None = None
    clade: str = ""


@dataclass
class GeneBox:
    gene_id: str
    name: str
    start: int
    end: int
    strand: str
    is_goi: bool


def _parse_gff_attr(field: str) -> dict[str, str]:
    out: dict[str, str] = {}
    for part in field.split(";"):
        if "=" in part:
            key, val = part.split("=", 1)
            out[key] = val
    return out


def _display_name(attrs: dict[str, str], gene_id: str) -> str:
    for key in ("Name", "gene", "ID"):
        if key in attrs and attrs[key]:
            text = attrs[key]
            if key == "ID" and text.startswith("gene-"):
                text = text.removeprefix("gene-")
            return text.replace("%2C", ",")[:28]
    return gene_id.removeprefix("gene-")[:28]


def parse_genes_on_chrom(gff_path: Path, chrom: str) -> list[GeneBox]:
    genes: list[GeneBox] = []
    if not gff_path.exists():
        return genes
    with gff_path.open(encoding="utf-8", errors="replace") as handle:
        for line in handle:
            if line.startswith("#") or "\tgene\t" not in line:
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9 or parts[0] != chrom:
                continue
            start, end = int(parts[3]), int(parts[4])
            strand = parts[6]
            attrs = _parse_gff_attr(parts[8])
            gid = attrs.get("ID", "")
            if not gid:
                continue
            name = _display_name(attrs, gid)
            genes.append(
                GeneBox(
                    gene_id=gid,
                    name=name,
                    start=start,
                    end=end,
                    strand=strand,
                    is_goi=False,
                )
            )
    genes.sort(key=lambda g: g.start)
    return genes


def _goi_span(genes: list[GeneBox], goi_ids: set[str]) -> tuple[int, int] | None:
    hits = [g for g in genes if g.gene_id in goi_ids or g.name in goi_ids]
    if not hits:
        return None
    return min(g.start for g in hits), max(g.end for g in hits)


def genes_in_window(
    all_genes: list[GeneBox],
    goi_ids: set[str],
    *,
    flank_bp: int = 180_000,
    max_genes: int = 22,
) -> list[GeneBox]:
    span = _goi_span(all_genes, goi_ids)
    if span is None:
        return []
    gstart, gend = span
    win_start = gstart - flank_bp
    win_end = gend + flank_bp
    selected = [g for g in all_genes if g.end >= win_start and g.start <= win_end]
    if len(selected) > max_genes:
        # Keep all GOI plus closest neighbours
        goi_list = [g for g in selected if g.gene_id in goi_ids or g.name in goi_ids]
        others = [g for g in selected if g not in goi_list]
        others.sort(key=lambda g: min(abs(g.start - gstart), abs(g.end - gend)))
        selected = sorted(goi_list + others[: max_genes - len(goi_list)], key=lambda g: g.start)
    return [
        GeneBox(
            gene_id=g.gene_id,
            name=g.name,
            start=g.start,
            end=g.end,
            strand=g.strand,
            is_goi=g.gene_id in goi_ids or g.name in goi_ids,
        )
        for g in selected
    ]


def _innexin_gff(organism_folder: str) -> Path | None:
    root = ANNOTATIONS_DIR / "innexin" / organism_folder.replace(" ", "_")
    if not root.is_dir():
        return None
    matches = sorted(root.glob("*.gff"))
    return matches[0] if matches else None


def default_panels() -> list[tuple[str, str, list[LocusSpec]]]:
    """Predefined panels mirroring thesis case studies."""
    dmel_gff = _innexin_gff("Drosophila_melanogaster")
    aedes_gff = _innexin_gff("Aedes_aegypti")
    anopheles_gff = _innexin_gff("Anopheles_gambiae")
    cele_gff = _innexin_gff("Caenorhabditis_elegans")
    schisto_gff = _innexin_gff("Schistocerca_americana")

    panels: list[tuple[str, str, list[LocusSpec]]] = []

    if dmel_gff:
        panels.append(
            (
                "dmel_innexin_cluster_proximal",
                "D. melanogaster innexin cluster — proximal locus (~6.9 Mb)",
                [
                    LocusSpec(
                        "D. melanogaster — ogre",
                        dmel_gff,
                        "NC_004354.4",
                        ("gene-Dmel_CG3039", "ogre"),
                        6_840_000,
                        7_020_000,
                        clade="insect",
                    ),
                    LocusSpec(
                        "D. melanogaster — Inx2",
                        dmel_gff,
                        "NC_004354.4",
                        ("gene-Dmel_CG4590", "Inx2"),
                        6_840_000,
                        7_020_000,
                        clade="insect",
                    ),
                ],
            )
        )
        panels.append(
            (
                "dmel_innexin_cluster_distal",
                "D. melanogaster innexin cluster — distal locus (~20.7 Mb)",
                [
                    LocusSpec(
                        "D. melanogaster — Inx2 (distal)",
                        dmel_gff,
                        "NC_004354.4",
                        ("gene-Dmel_CG4590", "Inx2"),
                        20_720_000,
                        20_795_000,
                        clade="insect",
                    ),
                    LocusSpec(
                        "D. melanogaster — Inx3",
                        dmel_gff,
                        "NT_033777.3",
                        ("gene-Dmel_CG1448", "Inx3"),
                        clade="insect",
                    ),
                ],
            )
        )

    cross_rows: list[LocusSpec] = []
    if dmel_gff:
        cross_rows.append(
            LocusSpec(
                "D. melanogaster — Inx2",
                dmel_gff,
                "NC_004354.4",
                ("gene-Dmel_CG4590", "Inx2"),
                clade="insect",
            )
        )
    if aedes_gff:
        cross_rows.append(
            LocusSpec(
                "Aedes aegypti — shakB",
                aedes_gff,
                "NC_035107.1",
                ("gene-LOC5566218", "shakB"),
                clade="insect",
            )
        )
    if anopheles_gff:
        cross_rows.append(
            LocusSpec(
                "Anopheles gambiae — shakB",
                anopheles_gff,
                "NC_064600.1",
                ("gene-LOC1272246", "shakB"),
                clade="insect",
            )
        )
    if schisto_gff:
        cross_rows.append(
            LocusSpec(
                "Schistocerca americana — inx2",
                schisto_gff,
                "NC_060120.1",
                ("gene-LOC124594079", "inx2"),
                clade="insect",
            )
        )
    if cele_gff:
        cross_rows.append(
            LocusSpec(
                "C. elegans — inx-2",
                cele_gff,
                "NC_003284.9",
                ("gene-CELE_F08G12.10", "inx-2"),
                clade="nematode",
            )
        )
    if cross_rows:
        panels.append(
            (
                "inx2_cross_clade_gene_order",
                "Inx2 / shakB / inx-2 locus — gene order by species",
                cross_rows,
            )
        )

    return panels


def load_row(spec: LocusSpec) -> tuple[LocusSpec, list[GeneBox]]:
    genes = parse_genes_on_chrom(spec.gff_path, spec.chrom)
    goi_ids = set(spec.goi_gene_ids)
    if spec.window_start is not None and spec.window_end is not None:
        windowed = [
            g for g in genes if g.end >= spec.window_start and g.start <= spec.window_end
        ]
        genes_out = [
            GeneBox(
                gene_id=g.gene_id,
                name=g.name,
                start=g.start,
                end=g.end,
                strand=g.strand,
                is_goi=g.gene_id in goi_ids or g.name in goi_ids,
            )
            for g in windowed
        ]
    else:
        genes_out = genes_in_window(genes, goi_ids)
    return spec, genes_out


def _wrap_label(text: str, width: int = 11) -> str:
    words = re.split(r"[\s\-]+", text)
    lines: list[str] = []
    current: list[str] = []
    for word in words:
        candidate = " ".join(current + [word])
        if len(candidate) <= width:
            current.append(word)
        else:
            if current:
                lines.append(" ".join(current))
            current = [word]
    if current:
        lines.append(" ".join(current))
    return "\n".join(lines[:3])


def _gene_color_key(name: str) -> str:
    return re.sub(r"^gene-", "", name.strip(), flags=re.IGNORECASE).casefold()


def _build_gene_color_map(rows: list[tuple[LocusSpec, list[GeneBox]]]) -> dict[str, str]:
    """Assign one color per gene name (case-insensitive); GOI kept red separately."""
    counts: dict[str, int] = {}
    for _, genes in rows:
        for gene in genes:
            if gene.is_goi:
                continue
            key = _gene_color_key(gene.name)
            counts[key] = counts.get(key, 0) + 1
    # Shared genes first so they get distinct, stable palette slots.
    ordered = sorted(counts, key=lambda k: (-counts[k], k))
    return {name: GENE_PALETTE[i % len(GENE_PALETTE)] for i, name in enumerate(ordered)}


def plot_panel(
    title: str,
    rows: list[tuple[LocusSpec, list[GeneBox]]],
    out_path: Path,
) -> None:
    if not rows or all(not genes for _, genes in rows):
        raise ValueError(f"No genes to plot for {title}")

    color_map = _build_gene_color_map(rows)
    n = len(rows)
    fig_h = max(2.8, 1.15 * n + 1.4)
    fig, axes = plt.subplots(n, 1, figsize=(14, fig_h), squeeze=False)

    for ax, (spec, genes) in zip(axes.flat, rows):
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis("off")

        if not genes:
            ax.text(0.5, 0.5, "No genes in window", ha="center", va="center")
            continue

        n_genes = len(genes)
        gap = 0.012
        box_w = (1.0 - ROW_LABEL_W - gap * (n_genes - 1)) / n_genes
        box_w = min(box_w, 0.055)
        total_w = n_genes * box_w + (n_genes - 1) * gap
        x0 = ROW_LABEL_W + (1.0 - ROW_LABEL_W - total_w) / 2

        strand = genes[0].strand if genes else "?"
        meta = f"{spec.chrom}\n{strand} strand"
        ax.text(
            0.01,
            0.5,
            f"{spec.species_label}\n{meta}",
            va="center",
            ha="left",
            fontsize=8,
            fontweight="bold",
            linespacing=1.25,
        )

        for i, gene in enumerate(genes):
            x = x0 + i * (box_w + gap)
            if gene.is_goi:
                color = GOI_COLOR
                edge = "#7B241C"
            else:
                color = color_map.get(_gene_color_key(gene.name), FLANK_COLOR)
                edge = "#2C3E50"
            rect = FancyBboxPatch(
                (x, 0.22),
                box_w,
                0.56,
                boxstyle="round,pad=0.012,rounding_size=0.02",
                linewidth=1.4 if gene.is_goi else 0.6,
                edgecolor=edge,
                facecolor=color,
                transform=ax.transAxes,
            )
            ax.add_patch(rect)
            label = _wrap_label(gene.name)
            ax.text(
                x + box_w / 2,
                0.5,
                label,
                ha="center",
                va="center",
                fontsize=6 if gene.is_goi else 5.5,
                fontweight="bold" if gene.is_goi else "normal",
                color="white" if gene.is_goi else "#1a1d26",
                transform=ax.transAxes,
                linespacing=0.95,
            )

    fig.suptitle(title, fontsize=13, fontweight="bold", y=0.98)
    fig.text(
        0.5,
        0.01,
        "Gene order along chromosome (NCBI GFF) · red = gene of interest · same name = same color",
        ha="center",
        fontsize=8,
        color="#555",
    )
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout(rect=(0, 0.03, 1, 0.95))
    fig.savefig(out_path, dpi=200, bbox_inches="tight")
    plt.close(fig)


def build_panel(panel_id: str, title: str, specs: list[LocusSpec], out_dir: Path) -> dict[str, str]:
    rows = [load_row(spec) for spec in specs]
    png = out_dir / f"{panel_id}.png"
    html = out_dir / f"{panel_id}.html"
    plot_panel(title, rows, png)

    # Simple HTML wrapper for gallery embedding
    rel_png = png.name
    html.write_text(
        f"""<!DOCTYPE html><html><head><meta charset='utf-8'>
<title>{title}</title>
<style>body{{font-family:system-ui,sans-serif;margin:2rem;background:#fafafa}}
img{{max-width:100%;border-radius:8px;box-shadow:0 4px 20px rgba(0,0,0,.1)}}
h1{{color:#2E86AB;font-size:1.2rem}}</style></head><body>
<h1>{title}</h1>
<p>GFF microsynteny — annotated gene order around gap-junction locus (no SynVoy).</p>
<img src='{rel_png}' alt='{title}'></body></html>""",
        encoding="utf-8",
    )

    manifest_rows = []
    for spec, genes in rows:
        for g in genes:
            manifest_rows.append(
                {
                    "panel": panel_id,
                    "species": spec.species_label,
                    "gene_id": g.gene_id,
                    "name": g.name,
                    "start": str(g.start),
                    "end": str(g.end),
                    "is_goi": "yes" if g.is_goi else "no",
                }
            )
    return {"panel": panel_id, "png": str(png), "html": str(html), "genes": str(len(manifest_rows))}


def build_all(out_dir: Path) -> list[dict[str, str]]:
    results: list[dict[str, str]] = []
    all_manifest: list[dict[str, str]] = []
    for panel_id, title, specs in default_panels():
        try:
            meta = build_panel(panel_id, title, specs, out_dir)
            results.append(meta)
            print(f"OK  {panel_id} -> {meta['png']}")
        except (ValueError, OSError) as exc:
            print(f"SKIP {panel_id}: {exc}", file=sys.stderr)
    manifest_path = out_dir / "gene_order_manifest.csv"
    for panel_id, _title, specs in default_panels():
        for spec in specs:
            _, genes = load_row(spec)
            for g in genes:
                all_manifest.append(
                    {
                        "panel": panel_id,
                        "species": spec.species_label,
                        "clade": spec.clade,
                        "gene_id": g.gene_id,
                        "name": g.name,
                        "chrom": spec.chrom,
                        "start": str(g.start),
                        "end": str(g.end),
                        "strand": g.strand,
                        "is_goi": "yes" if g.is_goi else "no",
                    }
                )
    write_csv(manifest_path, all_manifest)
    index = out_dir / "index.html"
    cards = "\n".join(
        f"<div class='card'><a href='{Path(r['html']).name}'><img src='{Path(r['png']).name}' "
        f"alt='{r['panel']}'><p>{r['panel'].replace('_', ' ')}</p></a></div>"
        for r in results
    )
    index.write_text(
        f"""<!DOCTYPE html><html><head><meta charset='utf-8'>
<title>Gene-order synteny</title>
<style>
body{{font-family:system-ui,sans-serif;max-width:1100px;margin:auto;padding:2rem;background:#f4f7fb}}
h1{{color:#2E86AB}}.grid{{display:grid;grid-template-columns:1fr;gap:2rem}}
.card img{{width:100%;border-radius:8px}}.card p{{text-align:center;font-weight:600}}
a{{text-decoration:none;color:inherit}}
</style></head><body>
<h1>Gene-order microsynteny (innexin loci)</h1>
<p>Full gene-order context from NCBI annotations (±180 kb neighbour windows).</p>
<div class='grid'>{cards}</div></body></html>""",
        encoding="utf-8",
    )
    return results


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build gene-order synteny plots from GFF")
    parser.add_argument("--outdir", type=Path, default=OUT_DIR)
    parser.add_argument("--panel", action="append", default=[], help="Panel id (default: all)")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> int:
    args = parse_args(argv)
    args.outdir.mkdir(parents=True, exist_ok=True)
    if args.panel:
        panels = default_panels()
        wanted = {p.lower() for p in args.panel}
        panels = [p for p in panels if p[0].lower() in wanted]
        if not panels:
            print("No panels matched.", file=sys.stderr)
            return 1
        for panel_id, title, specs in panels:
            try:
                build_panel(panel_id, title, specs, args.outdir)
                print(f"OK  {panel_id}")
            except (ValueError, OSError) as exc:
                print(f"FAIL {panel_id}: {exc}", file=sys.stderr)
                return 1
        return 0
    build_all(args.outdir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
