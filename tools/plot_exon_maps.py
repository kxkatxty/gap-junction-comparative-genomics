from __future__ import annotations

import gzip
import re
import subprocess
from pathlib import Path

import matplotlib.pyplot as plt
import pandas as pd
from matplotlib.patches import Rectangle
from matplotlib.transforms import blended_transform_factory

from fetch_protein_features import (
    FEATURES_CSV,
    WHITELIST_CSV,
    load_mapping,
    load_whitelist,
    match_whitelist_entries,
)


PROJECT_ROOT = Path("project")
RESULTS_DIR = PROJECT_ROOT / "results" / "exon_structures"
BIOCENTRAL_CSV = RESULTS_DIR / "protein_structure_predictions.csv"
ANNOTATION_DIR = PROJECT_ROOT / "data" / "annotations"
PLOTS_DIR = RESULTS_DIR / "plots"
FEATURE_PLOTS_DIR = PLOTS_DIR / "protein_feature_maps"
PLOTS_DIR.mkdir(parents=True, exist_ok=True)
FEATURE_PLOTS_DIR.mkdir(parents=True, exist_ok=True)

MAX_ROWS_PER_FIGURE = 80
ROW_HEIGHT = 0.22
TRACK_BAR_HEIGHT = 0.24
PROTEIN_BAR_Y = 0.03
CLASSIC_TRACK_CENTERS = {
    "pfam": 0.88,
    "coiled_coil": 0.64,
    "disorder": 0.40,
    "exons": 0.16,
}
ROSTLAB_TRACK_CENTERS = {
    "pfam": 0.90,
    "secondary_structure": 0.75,
    "transmembrane": 0.60,
    "disorder": 0.45,
    "coiled_coil": 0.30,
    "exons": 0.14,
}
FEATURE_FIG_WIDTH = 22
FEATURE_PANEL_HEIGHT = 3.8
FEATURE_DPI = 250
FEATURE_LEFT_MARGIN = 0.15
TRACK_LABEL_X = -0.018
SPECIES_LABEL_FONTSIZE = 18
SPECIES_FIG_X = 0.02
TRACK_LABEL_FONTSIZE = 13
DOMAIN_LABEL_FONTSIZE = 13
BRIGHTER_COLORS = {
    "#E67E22": "#FFB020",
    "#8E44AD": "#D68CFF",
    "#2980B9": "#66C2FF",
    "#C0392B": "#FF7070",
    "#27AE60": "#44E68A",
    "#BDC3C7": "#F2F2F2",
}
CLASSIC_TRACK_ORDER = ("pfam", "coiled_coil", "disorder", "exons")
ROSTLAB_TRACK_ORDER = (
    "pfam",
    "secondary_structure",
    "transmembrane",
    "disorder",
    "coiled_coil",
    "exons",
)
TRACK_LABELS = {
    "pfam": "Pfam",
    "secondary_structure": "Sec. struct.",
    "transmembrane": "TM topology",
    "coiled_coil": "Coiled",
    "disorder": "Disorder",
    "exons": "Exons",
}

_TRANSCRIPT_FEATURE_CACHE: dict[tuple[str, str], tuple[list[tuple[int, int]], list[tuple[int, int, int]]]] = {}
_LARGE_GFF_BYTES = 50_000_000


def normalize_species_name(name: str) -> str:
    return name.strip().replace(" ", "_")


def parse_attributes(attr_str: str) -> dict[str, str]:
    attrs: dict[str, str] = {}
    if "=" in attr_str:
        for field in attr_str.split(";"):
            field = field.strip()
            if not field or "=" not in field:
                continue
            key, value = field.split("=", 1)
            attrs[key.strip()] = value.strip().strip('"')
    else:
        for match in re.finditer(r'(\S+)\s+"([^"]+)"', attr_str):
            attrs[match.group(1)] = match.group(2)
    return attrs


def find_annotation_file(family: str, organism: str) -> Path | None:
    base = ANNOTATION_DIR / family.lower() / normalize_species_name(organism)
    if not base.exists():
        return None

    candidates: list[Path] = []
    for pattern in ("*.gff", "*.gff3", "*.gtf", "*.gff.gz", "*.gff3.gz"):
        candidates.extend(base.glob(pattern))
    if not candidates:
        return None

    def sort_key(path: Path) -> tuple[int, str]:
        name = path.name.lower()
        if name.endswith(".gff3.gz"):
            return 0, path.name
        if name.endswith(".gff3"):
            return 1, path.name
        if name.endswith(".gff.gz"):
            return 2, path.name
        if name.endswith(".gff"):
            return 3, path.name
        return 4, path.name

    return sorted(candidates, key=sort_key)[0]


def open_annotation(path: Path):
    if path.suffix.lower() == ".gz" or path.name.lower().endswith(".gff.gz"):
        return gzip.open(path, "rt", encoding="utf-8")
    return path.open("r", encoding="utf-8")


def _parse_feature_line(
    line: str, transcript_id: str
) -> tuple[str, tuple[int, int]] | tuple[str, tuple[int, int, int]] | None:
    if line.startswith("#"):
        return None
    parts = line.rstrip("\n").split("\t")
    if len(parts) != 9:
        return None

    _seqid, _source, feature_type, start, end, _score, _strand, phase, attributes = parts
    attrs = parse_attributes(attributes)
    parents = attrs.get("Parent") or attrs.get("transcript_id") or ""
    parent_ids = [p for p in re.split(r"[,\s]+", parents.strip()) if p]
    if transcript_id not in parent_ids:
        return None

    start_i = int(float(start))
    end_i = int(float(end))
    feature = feature_type.lower()
    if feature == "exon":
        return ("exon", (start_i, end_i))
    if feature == "cds":
        phase_i = 0 if phase == "." else int(phase)
        return ("cds", (start_i, end_i, phase_i))
    return None


def _grep_transcript_features(
    annotation_path: Path, transcript_id: str
) -> tuple[list[tuple[int, int]], list[tuple[int, int, int]]]:
    needle = f"Parent={transcript_id}"
    if annotation_path.name.lower().endswith(".gz"):
        command = ["zgrep", "-F", needle, str(annotation_path)]
    else:
        command = ["grep", "-F", needle, str(annotation_path)]

    exons: list[tuple[int, int]] = []
    cds_segments: list[tuple[int, int, int]] = []
    result = subprocess.run(command, capture_output=True, text=True, check=False)
    for line in result.stdout.splitlines():
        parsed = _parse_feature_line(line, transcript_id)
        if not parsed:
            continue
        kind, payload = parsed
        if kind == "exon":
            exons.append(payload)
        else:
            cds_segments.append(payload)
    return exons, cds_segments


def parse_transcript_cds_and_exons(
    annotation_path: Path, transcript_id: str
) -> tuple[list[tuple[int, int]], list[tuple[int, int, int]]]:
    cache_key = (str(annotation_path.resolve()), transcript_id)
    if cache_key in _TRANSCRIPT_FEATURE_CACHE:
        return _TRANSCRIPT_FEATURE_CACHE[cache_key]

    if annotation_path.stat().st_size >= _LARGE_GFF_BYTES:
        parsed = _grep_transcript_features(annotation_path, transcript_id)
    else:
        exons: list[tuple[int, int]] = []
        cds_segments: list[tuple[int, int, int]] = []
        with open_annotation(annotation_path) as handle:
            for line in handle:
                parsed = _parse_feature_line(line, transcript_id)
                if not parsed:
                    continue
                kind, payload = parsed
                if kind == "exon":
                    exons.append(payload)
                else:
                    cds_segments.append(payload)
        parsed = (exons, cds_segments)

    _TRANSCRIPT_FEATURE_CACHE[cache_key] = parsed
    return parsed


def read_fasta_sequence(path: Path) -> str:
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                chunks.append(line.strip())
    return "".join(chunks)


def map_cds_to_aa(
    cds_segments: list[tuple[int, int, int]], strand: str
) -> tuple[list[tuple[int, int, int, int]], int]:
    if not cds_segments:
        return [], 0

    ordered = (
        sorted(cds_segments, key=lambda item: -item[0])
        if strand == "-"
        else sorted(cds_segments, key=lambda item: item[0])
    )

    aa = 1
    pieces: list[tuple[int, int, int, int]] = []
    for index, (start, end, phase) in enumerate(ordered):
        length = end - start + 1
        if index == 0 and phase:
            usable = length - phase
            if usable <= 0:
                continue
            aa_len = usable // 3
            if aa_len <= 0:
                continue
            if strand == "+":
                coding_start = start + phase
                coding_end = coding_start + aa_len * 3 - 1
            else:
                coding_end = end - phase
                coding_start = coding_end - aa_len * 3 + 1
        else:
            aa_len = length // 3
            if aa_len <= 0:
                continue
            if strand == "+":
                coding_start = start
                coding_end = start + aa_len * 3 - 1
            else:
                coding_end = end
                coding_start = end - aa_len * 3 + 1

        pieces.append((coding_start, coding_end, aa, aa + aa_len - 1))
        aa += aa_len

    return pieces, max(aa - 1, 0)


def genomic_overlap_to_aa(
    overlap_start: int,
    overlap_end: int,
    coding_start: int,
    aa_start: int,
) -> tuple[int, int]:
    offset_start = overlap_start - coding_start
    offset_end = overlap_end - coding_start
    return aa_start + offset_start // 3, aa_start + offset_end // 3


def coding_exon_aa_ranges(
    exons: list[tuple[int, int]],
    cds_segments: list[tuple[int, int, int]],
    strand: str,
    protein_length: int,
) -> list[dict]:
    pieces, cds_aa_length = map_cds_to_aa(cds_segments, strand)
    if not pieces:
        return []

    ordered_exons = (
        sorted(exons, key=lambda item: -item[0])
        if strand == "-"
        else sorted(exons, key=lambda item: item[0])
    )

    coding_exons: list[dict] = []
    coding_number = 0
    for exon_start, exon_end in ordered_exons:
        aa_ranges: list[tuple[int, int]] = []
        for coding_start, coding_end, aa_start, aa_end in pieces:
            overlap_start = max(exon_start, coding_start)
            overlap_end = min(exon_end, coding_end)
            if overlap_start > overlap_end:
                continue
            ov_aa_start, ov_aa_end = genomic_overlap_to_aa(
                overlap_start, overlap_end, coding_start, aa_start
            )
            aa_ranges.append((ov_aa_start, ov_aa_end))

        if not aa_ranges:
            continue

        coding_number += 1
        coding_exons.append(
            {
                "coding_exon_number": coding_number,
                "aa_start": min(item[0] for item in aa_ranges),
                "aa_end": max(item[1] for item in aa_ranges),
            }
        )

    if protein_length > 0 and cds_aa_length > 0 and cds_aa_length != protein_length:
        scale = protein_length / cds_aa_length
        for exon in coding_exons:
            exon["aa_start"] = max(1, int(round((exon["aa_start"] - 1) * scale + 1)))
            exon["aa_end"] = max(exon["aa_start"], int(round(exon["aa_end"] * scale)))

    return coding_exons


def label_for_row(row: pd.Series) -> str:
    ref_gene = str(row.get("reference_gene_symbol", "") or "").strip()
    gff_gene = str(row.get("gene_symbol", "") or "").strip()
    organism = str(row.get("organism", "") or "").strip()

    if ref_gene and ref_gene != gff_gene:
        return f"{organism} | {ref_gene} ({gff_gene})"
    if ref_gene:
        return f"{organism} | {ref_gene}"
    return f"{organism} | {gff_gene}"


def pick_representative_transcripts(exon_df: pd.DataFrame) -> pd.DataFrame:
    if exon_df.empty:
        return exon_df

    transcript_stats = (
        exon_df.groupby(
            ["family", "organism", "fasta_path", "uniprot_accession", "transcript_id"],
            dropna=False,
        )
        .agg(
            reference_gene_symbol=("reference_gene_symbol", "first"),
            gene_symbol=("gene_symbol", "first"),
            transcript_span=("exon_end", lambda s: s.max() - s.min()),
            exon_count=("exon_number", "count"),
        )
        .reset_index()
    )

    transcript_stats = transcript_stats.sort_values(
        ["family", "organism", "fasta_path", "exon_count", "transcript_span"],
        ascending=[True, True, True, False, False],
    )
    chosen = transcript_stats.drop_duplicates(
        subset=["family", "fasta_path"], keep="first"
    )

    return exon_df.merge(
        chosen[["family", "fasta_path", "transcript_id"]],
        on=["family", "fasta_path", "transcript_id"],
        how="inner",
    )


def transcript_keys(df: pd.DataFrame) -> list[tuple]:
    keys = (
        df[
            [
                "organism",
                "fasta_path",
                "transcript_id",
                "reference_gene_symbol",
                "gene_symbol",
            ]
        ]
        .drop_duplicates()
        .sort_values(["organism", "reference_gene_symbol", "fasta_path"])
    )
    return list(keys.itertuples(index=False, name=None))


def chunk_keys(keys: list[tuple], chunk_size: int) -> list[list[tuple]]:
    return [keys[i : i + chunk_size] for i in range(0, len(keys), chunk_size)]


def get_coding_exons_for_transcript(
    row: pd.Series, tsub: pd.DataFrame
) -> tuple[list[dict], int]:
    fasta_path = Path(str(row["fasta_path"]))
    protein_length = len(read_fasta_sequence(fasta_path)) if fasta_path.exists() else 0

    annotation_path = find_annotation_file(str(row["family"]), str(row["organism"]))
    if not annotation_path:
        return [], protein_length

    transcript_id = str(tsub["transcript_id"].iloc[0])
    strand = str(tsub["strand"].iloc[0])
    gff_exons, cds_segments = parse_transcript_cds_and_exons(
        annotation_path, transcript_id
    )
    if not gff_exons:
        gff_exons = [
            (int(r["exon_start"]), int(r["exon_end"]))
            for _, r in tsub.sort_values("exon_number").iterrows()
        ]

    coding_exons = coding_exon_aa_ranges(
        gff_exons, cds_segments, strand, protein_length
    )
    return coding_exons, protein_length


def plot_compressed_exon_map(df: pd.DataFrame, family: str, out_file: Path) -> None:
    sub = pick_representative_transcripts(df[df["family"] == family].copy())
    if sub.empty:
        return

    keys = transcript_keys(sub)
    for part_index, key_chunk in enumerate(chunk_keys(keys, MAX_ROWS_PER_FIGURE), start=1):
        suffix = f"_part{part_index}" if len(keys) > MAX_ROWS_PER_FIGURE else ""
        part_file = out_file.with_name(out_file.stem + suffix + out_file.suffix)

        fig_h = max(4, min(len(key_chunk) * ROW_HEIGHT + 1.5, 36))
        fig, ax = plt.subplots(figsize=(14, fig_h))

        y = 0
        for organism, fasta_path, transcript_id, reference_gene_symbol, gene_symbol in key_chunk:
            tsub = sub[
                (sub["organism"] == organism)
                & (sub["fasta_path"] == fasta_path)
                & (sub["transcript_id"] == transcript_id)
            ].sort_values("exon_number")

            x = 0
            for _, erow in tsub.iterrows():
                width = max(0.4, erow["exon_length"] / 200)
                rect = Rectangle(
                    (x, y - 0.2), width, 0.4, facecolor="lightgray", edgecolor="black"
                )
                ax.add_patch(rect)
                ax.text(
                    x + width / 2,
                    y,
                    str(int(erow["exon_number"])),
                    ha="center",
                    va="center",
                    fontsize=7,
                )
                x += width + 0.25

            label_row = pd.Series(
                {
                    "organism": organism,
                    "reference_gene_symbol": reference_gene_symbol,
                    "gene_symbol": gene_symbol,
                }
            )
            ax.text(-0.3, y, label_for_row(label_row), ha="right", va="center", fontsize=8)
            y += 1

        title_suffix = f" ({part_index})" if len(keys) > MAX_ROWS_PER_FIGURE else ""
        ax.set_title(f"Compressed Exon Map: {family}{title_suffix}")
        ax.set_yticks([])
        ax.set_xlabel("Compressed exon structure (one transcript per FASTA)")
        ax.set_xlim(-1, None)
        ax.set_ylim(-1, y)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        plt.tight_layout()
        plt.savefig(part_file, dpi=200)
        plt.close()


def plot_genomic_exon_map(df: pd.DataFrame, family: str, out_file: Path) -> None:
    sub = pick_representative_transcripts(df[df["family"] == family].copy())
    if sub.empty:
        return

    keys = transcript_keys(sub)
    for part_index, key_chunk in enumerate(chunk_keys(keys, MAX_ROWS_PER_FIGURE), start=1):
        suffix = f"_part{part_index}" if len(keys) > MAX_ROWS_PER_FIGURE else ""
        part_file = out_file.with_name(out_file.stem + suffix + out_file.suffix)

        fig_h = max(4, min(len(key_chunk) * ROW_HEIGHT + 1.5, 36))
        fig, ax = plt.subplots(figsize=(16, fig_h))

        y = 0
        xmax = 0
        for organism, fasta_path, transcript_id, reference_gene_symbol, gene_symbol in key_chunk:
            tsub = sub[
                (sub["organism"] == organism)
                & (sub["fasta_path"] == fasta_path)
                & (sub["transcript_id"] == transcript_id)
            ].sort_values("exon_number")

            tx_start = tsub["exon_start"].min()
            tx_end = tsub["exon_end"].max()
            rel_end = tx_end - tx_start
            ax.plot([0, rel_end], [y, y], color="black", linewidth=0.8)

            for _, row in tsub.iterrows():
                x0 = row["exon_start"] - tx_start
                width = row["exon_end"] - row["exon_start"] + 1
                rect = Rectangle(
                    (x0, y - 0.2), width, 0.4, facecolor="gray", edgecolor="black"
                )
                ax.add_patch(rect)

            label_row = pd.Series(
                {
                    "organism": organism,
                    "reference_gene_symbol": reference_gene_symbol,
                    "gene_symbol": gene_symbol,
                }
            )
            ax.text(-rel_end * 0.02, y, label_for_row(label_row), ha="right", va="center", fontsize=8)
            xmax = max(xmax, rel_end)
            y += 1

        title_suffix = f" ({part_index})" if len(keys) > MAX_ROWS_PER_FIGURE else ""
        ax.set_title(f"Genomic Exon Map: {family}{title_suffix}")
        ax.set_yticks([])
        ax.set_xlabel("Relative genomic position within transcript (bp)")
        ax.set_ylim(-1, y)
        ax.set_xlim(-xmax * 0.08, xmax * 1.05 if xmax else 1)
        ax.spines["left"].set_visible(False)
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        plt.tight_layout()
        plt.savefig(part_file, dpi=200)
        plt.close()


def _brighten_color(color: str) -> str:
    return BRIGHTER_COLORS.get(str(color), color)


def short_organism(name: str) -> str:
    parts = name.split()
    if len(parts) >= 2:
        return f"{parts[0][0]}. {' '.join(parts[1:])}"
    return name


def _blended_transform(ax):
    return blended_transform_factory(ax.transData, ax.transAxes)


def combine_protein_features(
    features_df: pd.DataFrame,
    biocentral_df: pd.DataFrame,
    accession: str,
) -> tuple[pd.DataFrame, bool]:
    accession = accession.upper()
    classic = (
        features_df[features_df["uniprot_accession"].astype(str).str.upper() == accession]
        if not features_df.empty
        else pd.DataFrame()
    )
    rostlab = (
        biocentral_df[biocentral_df["uniprot_accession"].astype(str).str.upper() == accession]
        if not biocentral_df.empty
        else pd.DataFrame()
    )
    has_rostlab = not rostlab.empty

    if not has_rostlab:
        return classic, False

    if classic.empty:
        return rostlab, True

    classic = classic.copy()
    if "disorder" in rostlab["track"].values:
        classic = classic[classic["track"] != "disorder"]
    return pd.concat([classic, rostlab], ignore_index=True), True


def _draw_feature_boxes(
    ax,
    features: pd.DataFrame,
    track: str,
    track_centers: dict[str, float],
    label_key: str = "feature_label",
) -> None:
    track_features = features[features["track"] == track]
    if track_features.empty:
        return

    y_center = track_centers[track]
    trans = _blended_transform(ax)
    for _, feat in track_features.iterrows():
        width = feat["aa_end"] - feat["aa_start"] + 1
        rect = Rectangle(
            (feat["aa_start"], y_center - TRACK_BAR_HEIGHT / 2),
            width,
            TRACK_BAR_HEIGHT,
            transform=trans,
            facecolor=_brighten_color(feat.get("color", "#CCCCCC")),
            edgecolor="black",
            linewidth=2.0,
            zorder=3,
            clip_on=False,
        )
        ax.add_patch(rect)
        if width >= 14:
            ax.text(
                feat["aa_start"] + width / 2,
                y_center,
                str(feat[label_key]),
                transform=trans,
                ha="center",
                va="center",
                fontsize=DOMAIN_LABEL_FONTSIZE,
                fontweight="bold",
                color="white" if track in {"pfam", "secondary_structure", "transmembrane"} else "#111111",
                zorder=4,
                clip_on=False,
            )


def _draw_species_feature_panel(
    ax,
    protein_length: int,
    coding_exons: list[dict],
    features: pd.DataFrame,
    show_xlabel: bool,
    *,
    track_order: tuple[str, ...],
    track_centers: dict[str, float],
) -> None:
    trans = _blended_transform(ax)

    for exon in coding_exons:
        for boundary in (exon["aa_start"], exon["aa_end"] + 1):
            ax.axvline(
                boundary,
                color="#D0D0D0",
                linewidth=0.7,
                zorder=1,
            )

    for track in track_order:
        if track == "exons":
            continue
        _draw_feature_boxes(ax, features, track, track_centers)

    exon_y = track_centers["exons"]
    for exon in coding_exons:
        width = exon["aa_end"] - exon["aa_start"] + 1
        rect = Rectangle(
            (exon["aa_start"], exon_y - TRACK_BAR_HEIGHT / 2),
            width,
            TRACK_BAR_HEIGHT,
            transform=trans,
            facecolor="white",
            edgecolor="black",
            linewidth=1.2,
            zorder=3,
            clip_on=False,
        )
        ax.add_patch(rect)
        ax.text(
            exon["aa_start"] + width / 2,
            exon_y,
            str(exon["coding_exon_number"]),
            transform=trans,
            ha="center",
            va="center",
            fontsize=DOMAIN_LABEL_FONTSIZE,
            fontweight="bold",
            zorder=4,
            clip_on=False,
        )

    ax.plot(
        [1, protein_length],
        [PROTEIN_BAR_Y, PROTEIN_BAR_Y],
        transform=trans,
        color="black",
        linewidth=5.0,
        zorder=2,
        clip_on=False,
    )

    ax.set_xlim(0.5, max(protein_length, 10) + 12)
    ax.set_ylim(0, 1)
    ax.set_yticks([])
    ax.tick_params(axis="x", labelsize=12)
    ax.grid(axis="x", color="#EEEEEE", linewidth=0.6, zorder=0)
    ax.spines["left"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.spines["top"].set_visible(False)

    for track_name in track_order:
        ax.text(
            TRACK_LABEL_X,
            track_centers[track_name],
            TRACK_LABELS[track_name],
            transform=ax.transAxes,
            ha="right",
            va="center",
            fontsize=TRACK_LABEL_FONTSIZE,
            fontweight="bold",
            clip_on=False,
        )

    if show_xlabel:
        ax.set_xlabel("Protein position (aa)", fontsize=13, labelpad=8)
    else:
        ax.tick_params(axis="x", labelbottom=False)


def _add_species_labels(fig, plotted_axes: list[tuple]) -> None:
    for ax, organism in plotted_axes:
        pos = ax.get_position()
        fig.text(
            SPECIES_FIG_X,
            pos.y0 + pos.height / 2,
            short_organism(organism),
            ha="left",
            va="center",
            fontsize=SPECIES_LABEL_FONTSIZE,
            fontstyle="italic",
            fontweight="bold",
        )


def _render_feature_map_group(
    *,
    plot_group: str,
    title_label: str,
    group_entries: pd.DataFrame,
    exon_sub: pd.DataFrame,
    features_df: pd.DataFrame,
    biocentral_df: pd.DataFrame,
    use_rostlab: bool,
) -> Path | None:
    track_order = ROSTLAB_TRACK_ORDER if use_rostlab else CLASSIC_TRACK_ORDER
    track_centers = ROSTLAB_TRACK_CENTERS if use_rostlab else CLASSIC_TRACK_CENTERS
    panel_h = FEATURE_PANEL_HEIGHT + (0.8 if use_rostlab else 0.0)

    group_entries = group_entries.sort_values("organism")
    n_species = len(group_entries)
    fig_h = max(8.0, panel_h * n_species + 1.8)
    fig, axes = plt.subplots(
        n_species, 1, figsize=(FEATURE_FIG_WIDTH, fig_h), squeeze=False
    )
    if use_rostlab:
        subtitle = (
            f"{title_label}: Pfam + Biocentral structure predictions "
            "(ProtT5, TMbed, SETH) with coding exon boundaries"
        )
        out_name = f"{plot_group}_structure_map.png"
    else:
        subtitle = (
            f"{title_label}: domains/features with coding exon boundaries "
            "(InterPro/UniProt-derived)"
        )
        out_name = f"{plot_group}_protein_feature_map.png"

    fig.suptitle(
        subtitle,
        fontsize=14,
        x=FEATURE_LEFT_MARGIN,
        y=0.985,
        ha="left",
    )

    plotted = 0
    plotted_axes: list[tuple] = []
    entry_list = list(group_entries.iterrows())
    for panel_index, (ax, (_, entry)) in enumerate(zip(axes.ravel(), entry_list)):
        tsub = exon_sub[
            (exon_sub["fasta_path"] == entry["fasta_path"])
            & (exon_sub["organism"] == entry["organism"])
        ]
        if tsub.empty:
            ax.axis("off")
            ax.text(0.5, 0.5, f"No exon data for {entry['organism']}", ha="center")
            continue

        tsub = tsub[tsub["transcript_id"] == tsub["transcript_id"].iloc[0]]
        coding_exons, protein_length = get_coding_exons_for_transcript(entry, tsub)
        if not coding_exons or protein_length <= 0:
            ax.axis("off")
            ax.text(
                0.5,
                0.5,
                f"CDS mapping unavailable for {entry['organism']}",
                ha="center",
            )
            continue

        accession = str(entry["uniprot_accession"]).upper()
        if use_rostlab:
            protein_features, _ = combine_protein_features(
                features_df, biocentral_df, accession
            )
        else:
            protein_features = features_df[
                features_df["uniprot_accession"].astype(str).str.upper() == accession
            ]

        _draw_species_feature_panel(
            ax,
            protein_length=protein_length,
            coding_exons=coding_exons,
            features=protein_features,
            show_xlabel=(panel_index == len(entry_list) - 1),
            track_order=track_order,
            track_centers=track_centers,
        )
        plotted += 1
        plotted_axes.append((ax, str(entry["organism"])))

    if plotted == 0:
        plt.close(fig)
        return None

    fig.subplots_adjust(
        left=FEATURE_LEFT_MARGIN,
        right=0.99,
        top=0.97,
        bottom=0.05,
        hspace=0.42,
    )
    _add_species_labels(fig, plotted_axes)
    out_file = FEATURE_PLOTS_DIR / out_name
    fig.savefig(out_file, dpi=FEATURE_DPI, pad_inches=0.25)
    plt.close(fig)
    return out_file


def plot_protein_feature_maps(
    exon_df: pd.DataFrame,
    whitelist_entries: pd.DataFrame,
    features_df: pd.DataFrame,
    biocentral_df: pd.DataFrame | None = None,
) -> list[Path]:
    created: list[Path] = []
    if whitelist_entries.empty:
        return created

    biocentral_df = biocentral_df if biocentral_df is not None else pd.DataFrame()
    sub = pick_representative_transcripts(exon_df)

    for plot_group, group_entries in whitelist_entries.groupby("plot_group"):
        title_label = str(group_entries["title_label"].iloc[0])
        classic_file = _render_feature_map_group(
            plot_group=plot_group,
            title_label=title_label,
            group_entries=group_entries,
            exon_sub=sub,
            features_df=features_df,
            biocentral_df=biocentral_df,
            use_rostlab=False,
        )
        if classic_file is not None:
            created.append(classic_file)

        group_accessions = group_entries["uniprot_accession"].astype(str).str.upper()
        has_rostlab = (
            not biocentral_df.empty
            and biocentral_df["uniprot_accession"]
            .astype(str)
            .str.upper()
            .isin(group_accessions)
            .any()
        )
        if has_rostlab:
            rostlab_file = _render_feature_map_group(
                plot_group=plot_group,
                title_label=title_label,
                group_entries=group_entries,
                exon_sub=sub,
                features_df=features_df,
                biocentral_df=biocentral_df,
                use_rostlab=True,
            )
            if rostlab_file is not None:
                created.append(rostlab_file)

    return created


def plot_exon_count_boxplot(transcript_df: pd.DataFrame, out_file: Path) -> None:
    if transcript_df.empty:
        return

    plot_df = (
        transcript_df.sort_values("exon_count", ascending=False)
        .drop_duplicates(subset=["family", "fasta_path"], keep="first")
    )

    fig, ax = plt.subplots(figsize=(8, 6))
    plot_df.boxplot(column="exon_count", by="family", ax=ax)
    plt.suptitle("")
    ax.set_title("Exon count per reference protein by family")
    ax.set_ylabel("Exon count (representative transcript)")
    plt.tight_layout()
    plt.savefig(out_file, dpi=200)
    plt.close()


def main() -> None:
    exon_csv = RESULTS_DIR / "exon_coordinates.csv"
    transcript_csv = RESULTS_DIR / "transcript_summary.csv"

    if not exon_csv.exists():
        raise FileNotFoundError(f"Missing {exon_csv}")
    if not transcript_csv.exists():
        raise FileNotFoundError(f"Missing {transcript_csv}")

    exon_df = pd.read_csv(exon_csv)
    transcript_df = pd.read_csv(transcript_csv)

    required = {
        "family",
        "organism",
        "fasta_path",
        "gene_symbol",
        "transcript_id",
        "exon_number",
        "exon_start",
        "exon_end",
        "exon_length",
        "strand",
    }
    missing = required - set(exon_df.columns)
    if missing:
        raise ValueError(f"exon_coordinates.csv missing columns: {sorted(missing)}")

    families = sorted(exon_df["family"].dropna().unique().tolist())

    for family in families:
        plot_compressed_exon_map(
            exon_df,
            family,
            PLOTS_DIR / f"{family}_compressed_exon_map.png",
        )
        plot_genomic_exon_map(
            exon_df,
            family,
            PLOTS_DIR / f"{family}_genomic_exon_map.png",
        )

    plot_exon_count_boxplot(
        transcript_df,
        PLOTS_DIR / "exon_count_boxplot_by_family.png",
    )

    feature_plots: list[Path] = []
    if WHITELIST_CSV.exists() and (RESULTS_DIR / "protein_mapping.csv").exists():
        whitelist = load_whitelist(WHITELIST_CSV)
        mapping = load_mapping(RESULTS_DIR / "protein_mapping.csv")
        whitelist_entries = match_whitelist_entries(whitelist, mapping)
        if FEATURES_CSV.exists():
            features_df = pd.read_csv(FEATURES_CSV)
        else:
            features_df = pd.DataFrame()
        biocentral_df = (
            pd.read_csv(BIOCENTRAL_CSV) if BIOCENTRAL_CSV.exists() else pd.DataFrame()
        )
        feature_plots = plot_protein_feature_maps(
            exon_df, whitelist_entries, features_df, biocentral_df
        )

    created = sorted(PLOTS_DIR.glob("*.png")) + sorted(FEATURE_PLOTS_DIR.glob("*.png"))
    print(f"Plots saved in: {PLOTS_DIR}")
    for path in created:
        print(f"  - {path.relative_to(PROJECT_ROOT)}")
    if WHITELIST_CSV.exists() and not FEATURES_CSV.exists():
        print(
            f"\nTip: run fetch_protein_features.py to add Pfam/domain overlays "
            f"for whitelist proteins."
        )
    if WHITELIST_CSV.exists() and not BIOCENTRAL_CSV.exists():
        print(
            f"\nTip: run fetch_biocentral_predictions.py to add Rostlab/Biocentral "
            f"structure tracks (ProtT5, TMbed, SETH)."
        )


if __name__ == "__main__":
    main()
