#!/usr/bin/env python3
"""Clade × innexin-type comparison across reference DB + discovery + curator loci.

Builds diagrams and a website under project/results/innexin_clade_comparison/.
"""

from __future__ import annotations

import csv
import html
import re
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.patches import Patch

from pipeline.common import METADATA_DIR, PROJECT_ROOT, RESULTS_DIR, write_csv

OUT = RESULTS_DIR / "innexin_clade_comparison"
FIGS = OUT / "figures"
CAND = METADATA_DIR / "gap_junction_candidates_master.csv"
REF_PANEL = METADATA_DIR / "innexin_reference_panel.csv"
REF_ROOT = PROJECT_ROOT / "project/data/references/innexins"

# Visual palette (no purple-on-white / cream-terracotta defaults)
CLADE_COLORS = {
    "Insecta": "#0F766E",
    "Nematoda": "#B45309",
    "Rotifera": "#0369A1",
    "Mollusca": "#BE123C",
    "Arthropoda (non-insect)": "#15803D",
    "Annelida": "#A16207",
    "Ctenophora": "#0E7490",
    "Reference insect": "#115E59",
    "Reference nematode": "#9A3412",
    "Other": "#64748B",
}

TYPE_COLORS = {
    "Inx1/ogre": "#0F766E",
    "Inx2": "#0369A1",
    "Inx3": "#BE123C",
    "Inx4/zpg": "#15803D",
    "Inx5": "#A16207",
    "Inx6": "#0E7490",
    "Inx7": "#7C2D12",
    "shakB": "#1D4ED8",
    "unc-7": "#B45309",
    "unc-9": "#C2410C",
    "inx-2 (Cele)": "#92400E",
    "other nematode": "#78716C",
    "other/unknown": "#94A3B8",
}

SUBFAMILY_FROM_TYPE = {
    "shakB": "SF1_shakB",
    "Inx1/ogre": "SF2_Inx1_ogre",
    "Inx3": "SF3_Inx3_Inx7_nematode",
    "Inx7": "SF3_Inx3_Inx7_nematode",
    "unc-7": "SF3_Inx3_Inx7_nematode",
    "unc-9": "SF3_Inx3_Inx7_nematode",
    "inx-2 (Cele)": "SF3_Inx3_Inx7_nematode",
    "other nematode": "SF3_Inx3_Inx7_nematode",
    "Inx2": "SF4_Inx2_expansion",
    "Inx4/zpg": "SF4_Inx2_expansion",
    "Inx5": "SF4_Inx2_expansion",
    "Inx6": "SF4_Inx2_expansion",
    "other/unknown": "unassigned",
}

SPECIES_CLADE = {
    # insects (reference + known)
    "Drosophila_melanogaster": "Insecta",
    "Aedes_aegypti": "Insecta",
    "Anopheles_gambiae": "Insecta",
    "Schistocerca_americana": "Insecta",
    # nematodes
    "Caenorhabditis_elegans": "Nematoda",
    "Caenorhabditis_briggsae": "Nematoda",
    # rotifers
    "Adineta_vaga": "Rotifera",
    "Adineta_ricciae": "Rotifera",
    "Adineta_steineri": "Rotifera",
    "Brachionus_calyciflorus": "Rotifera",
    "Brachionus_koreanus": "Rotifera",
    "Brachionus_manjavacas": "Rotifera",
    "Rotaria_macrura": "Rotifera",
    # molluscs
    "Abra_alba": "Mollusca",
    "Abra_segmentum": "Mollusca",
    "Acanthocardia_echinata": "Mollusca",
    "Acanthochitona_discrepans": "Mollusca",
    # arthropods non-insect
    "Acartia_tonsa": "Arthropoda (non-insect)",
    "Aegaeobuthus_cyprius": "Arthropoda (non-insect)",
    "Aelurillus_cypriotus": "Arthropoda (non-insect)",
    "Agelena_orientalis": "Arthropoda (non-insect)",
    "Amaurobius_ferox": "Arthropoda (non-insect)",
    "Eriophyidae_sp.": "Arthropoda (non-insect)",
    "Triops_cancriformis": "Arthropoda (non-insect)",
    # annelid / ctenophore
    "Streblospio_benedicti": "Annelida",
    "Dryodora_glandiformis": "Ctenophora",
    "Mnemiopsis_leidyi": "Ctenophora",
}

ACCEPTED_RANKS = {
    "curator_present_innexin",
    "high_confidence",
    "accepted",
    "weak_manual_review",
    "probable",
}


plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.titlesize": 12,
        "axes.labelsize": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "font.family": "DejaVu Sans",
    }
)


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _f(x, default=0.0) -> float:
    try:
        if x is None or x == "":
            return default
        return float(x)
    except (TypeError, ValueError):
        return default


def _i(x, default=0) -> int:
    try:
        if x is None or x == "":
            return default
        return int(float(x))
    except (TypeError, ValueError):
        return default


def infer_clade(species_dir: str, organism: str = "") -> str:
    if species_dir in SPECIES_CLADE:
        return SPECIES_CLADE[species_dir]
    name = (organism or species_dir or "").lower().replace("_", " ")
    if any(k in name for k in ("drosophila", "aedes", "anopheles", "schistocerca", "apis", "tribolium")):
        return "Insecta"
    if any(k in name for k in ("caenorhabditis", "pristionchus", "brugia")):
        return "Nematoda"
    if any(k in name for k in ("adineta", "brachionus", "rotaria", "rotifer")):
        return "Rotifera"
    if any(k in name for k in ("abra", "acanthocardia", "acanthochitona", "mytilus", "octopus", "crassostrea")):
        return "Mollusca"
    if any(k in name for k in ("streblospio", "capitella", "platynereis", "annel")):
        return "Annelida"
    if any(k in name for k in ("dryodora", "mnemiopsis", "beroe", "cteno")):
        return "Ctenophora"
    if any(
        k in name
        for k in (
            "agelena",
            "amaurobius",
            "aelurillus",
            "aegaeobuthus",
            "acartia",
            "triops",
            "eriophy",
            "spider",
            "scorpion",
        )
    ):
        return "Arthropoda (non-insect)"
    return "Other"


def classify_reference_type(best_hit: str, gene_hint: str = "") -> str:
    text = f"{best_hit} {gene_hint}".upper()
    # Order matters for specificity
    rules = [
        (r"SHAKB|SHAKE.?B", "shakB"),
        (r"INX1_DROME|\bOGRE\b|INX1_SCHAM|INX-?1\b", "Inx1/ogre"),
        (r"INX7_DROME|INX-?7_DROME|\bINX7\b", "Inx7"),
        (r"INX6_DROME|\bINX6\b", "Inx6"),
        (r"INX5_DROME|\bINX5\b", "Inx5"),
        (r"INX4_DROME|\bZPG\b|INX4_ANOGA", "Inx4/zpg"),
        (r"INX3_DROME|\bINX3\b(?!_CAE)", "Inx3"),
        (r"INX2_DROME|INX2_SCHAM|\bINX2\b(?!_CAE)", "Inx2"),
        (r"UNC-?7|UNC7_CAEEL", "unc-7"),
        (r"UNC-?9|UNC9_CAEEL", "unc-9"),
        (r"INX2_CAEEL|INX-?2_CAEEL", "inx-2 (Cele)"),
        (r"INX\d+_CAE|EAT-?5|UNC-|CAEEL|CAEBR", "other nematode"),
    ]
    for pat, label in rules:
        if re.search(pat, text):
            return label
    if "INX" in text or "INNEXIN" in text:
        return "other/unknown"
    return "other/unknown"


def _is_kept(rank: str) -> bool:
    r = (rank or "").lower()
    if r == "rejected_false_positive":
        return False
    if r in {x.lower() for x in ACCEPTED_RANKS}:
        return True
    # keep medium/high non-rejected discovery hits
    return r not in ("", "rejected_false_positive") and "reject" not in r


def load_reference_rows() -> list[dict]:
    rows: list[dict] = []
    if not REF_ROOT.exists():
        return rows
    for fasta in sorted(REF_ROOT.rglob("*.fasta")):
        organism = fasta.parent.name.replace("_", " ")
        stem = fasta.stem  # e.g. Inx2__Q9V427 or ogre__P27716
        gene = stem.split("__")[0]
        acc = stem.split("__")[1] if "__" in stem else ""
        seq = []
        with fasta.open(encoding="utf-8", errors="replace") as fh:
            for line in fh:
                if line.startswith(">"):
                    continue
                seq.append(line.strip())
        aa = len("".join(seq))
        species_dir = fasta.parent.name
        clade = infer_clade(species_dir, organism)
        if clade == "Insecta":
            clade = "Reference insect"
        elif clade == "Nematoda":
            clade = "Reference nematode"
        rtype = classify_reference_type(stem, gene)
        rows.append(
            {
                "source": "reference_db",
                "species_dir": species_dir,
                "organism": organism,
                "clade": clade,
                "candidate_id": stem,
                "gene_label": gene,
                "reference_type": rtype,
                "subfamily": SUBFAMILY_FROM_TYPE.get(rtype, "unassigned"),
                "protein_length": aa,
                "exon_count": "",
                "tm_helix_count": "",
                "reference_identity": 1.0,
                "best_reference_hit": f"{gene}|{acc}" if acc else gene,
                "rank_category": "reference",
                "discovery_batch": "reference_db",
            }
        )
    return rows


def load_candidate_rows() -> list[dict]:
    rows: list[dict] = []
    if not CAND.exists():
        return rows
    with CAND.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("family") != "innexin":
                continue
            rank = r.get("rank_category") or ""
            if not _is_kept(rank):
                continue
            species_dir = r.get("species_dir") or ""
            organism = r.get("organism") or species_dir.replace("_", " ")
            hit = r.get("best_reference_hit") or ""
            rtype = classify_reference_type(hit, r.get("candidate_id") or "")
            batch = r.get("discovery_batch") or ""
            source = "curator_new" if batch == "gene_curator_probe" else "discovery_old"
            rows.append(
                {
                    "source": source,
                    "species_dir": species_dir,
                    "organism": organism,
                    "clade": infer_clade(species_dir, organism),
                    "candidate_id": r.get("candidate_id") or "",
                    "gene_label": "",
                    "reference_type": rtype,
                    "subfamily": SUBFAMILY_FROM_TYPE.get(rtype, "unassigned"),
                    "protein_length": _i(r.get("protein_length")),
                    "exon_count": _i(r.get("exon_count")),
                    "tm_helix_count": _i(r.get("tm_helix_count")),
                    "reference_identity": _f(r.get("reference_identity")),
                    "best_reference_hit": hit,
                    "rank_category": rank,
                    "discovery_batch": batch,
                }
            )
    return rows


def summarize(rows: list[dict]) -> dict:
    kept = [r for r in rows if r["source"] != "noise"]
    by_clade = Counter(r["clade"] for r in kept)
    by_type = Counter(r["reference_type"] for r in kept)
    by_source = Counter(r["source"] for r in kept)
    by_sf = Counter(r["subfamily"] for r in kept)
    species = {(r["clade"], r["species_dir"]) for r in kept}
    return {
        "n_loci": len(kept),
        "n_species": len({r["species_dir"] for r in kept}),
        "by_clade": by_clade,
        "by_type": by_type,
        "by_source": by_source,
        "by_sf": by_sf,
        "species_per_clade": Counter(c for c, _ in species),
    }


def plot_loci_by_clade(rows: list[dict], out: Path) -> None:
    sources = ["reference_db", "discovery_old", "curator_new"]
    labels = ["Reference DB", "Old discovery", "Curator new"]
    colors = ["#115E59", "#64748B", "#0EA5E9"]
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(clades))
    bottom = np.zeros(len(clades))
    for src, lab, col in zip(sources, labels, colors):
        vals = np.array([sum(1 for r in rows if r["clade"] == c and r["source"] == src) for c in clades], dtype=float)
        ax.bar(x, vals, bottom=bottom, label=lab, color=col, width=0.72)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=25, ha="right")
    ax.set_ylabel("Innexin loci")
    ax.set_title("Innexin loci by clade (reference + old + curator)")
    ax.legend(frameon=False)
    _save(fig, out)


def plot_heatmap_clade_type(rows: list[dict], out: Path) -> None:
    types = [t for t, _ in Counter(r["reference_type"] for r in rows).most_common()]
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    mat = np.zeros((len(clades), len(types)))
    idx_c = {c: i for i, c in enumerate(clades)}
    idx_t = {t: i for i, t in enumerate(types)}
    for r in rows:
        mat[idx_c[r["clade"]], idx_t[r["reference_type"]]] += 1
    fig, ax = plt.subplots(figsize=(max(10, 0.7 * len(types) + 4), max(5, 0.45 * len(clades) + 2)))
    im = ax.imshow(mat, aspect="auto", cmap="YlGnBu")
    ax.set_xticks(range(len(types)))
    ax.set_xticklabels(types, rotation=40, ha="right", fontsize=8)
    ax.set_yticks(range(len(clades)))
    ax.set_yticklabels(clades, fontsize=9)
    for i in range(mat.shape[0]):
        for j in range(mat.shape[1]):
            v = int(mat[i, j])
            if v:
                ax.text(j, i, str(v), ha="center", va="center", fontsize=7, color="#0f172a" if v < mat.max() * 0.6 else "white")
    ax.set_title("Clade × closest reference type (best-hit labels)")
    fig.colorbar(im, ax=ax, fraction=0.03, pad=0.02, label="Locus count")
    _save(fig, out)


def _boxplot_by(rows: list[dict], key: str, value: str, out: Path, title: str, ylabel: str, min_n: int = 3) -> None:
    groups: dict[str, list[float]] = defaultdict(list)
    for r in rows:
        v = r.get(value)
        if value == "reference_identity":
            if not v or float(v) <= 0:
                continue
            groups[str(r[key])].append(float(v) * (100 if float(v) <= 1.5 else 1))
        else:
            n = _f(v)
            if n > 0:
                groups[str(r[key])].append(n)
    labels = [k for k, vals in sorted(groups.items(), key=lambda kv: -len(kv[1])) if len(vals) >= min_n]
    if not labels:
        return
    data = [groups[k] for k in labels]
    fig, ax = plt.subplots(figsize=(max(8, 0.55 * len(labels) + 3), 5))
    bp = ax.boxplot(data, tick_labels=labels, patch_artist=True, showfliers=False)
    for patch, lab in zip(bp["boxes"], labels):
        color = CLADE_COLORS.get(lab) or TYPE_COLORS.get(lab) or "#64748B"
        patch.set_facecolor(color)
        patch.set_alpha(0.75)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=30, ha="right")
    _save(fig, out)


def plot_type_composition_per_clade(rows: list[dict], out: Path) -> None:
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    types = [t for t, _ in Counter(r["reference_type"] for r in rows).most_common()]
    fig, ax = plt.subplots(figsize=(12, 5.8))
    x = np.arange(len(clades))
    bottom = np.zeros(len(clades))
    for t in types:
        vals = np.array([sum(1 for r in rows if r["clade"] == c and r["reference_type"] == t) for c in clades], dtype=float)
        if vals.sum() == 0:
            continue
        ax.bar(x, vals, bottom=bottom, label=t, color=TYPE_COLORS.get(t, "#94A3B8"), width=0.75)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=25, ha="right")
    ax.set_ylabel("Loci")
    ax.set_title("Reference-type mix within each clade")
    ax.legend(ncol=3, fontsize=7, frameon=False, loc="upper right")
    _save(fig, out)


def plot_species_richness(rows: list[dict], out: Path) -> None:
    # loci per species, colored by clade
    counts = Counter((r["clade"], r["organism"] or r["species_dir"]) for r in rows if r["source"] != "reference_db")
    items = sorted(counts.items(), key=lambda kv: -kv[1])[:30]
    if not items:
        return
    fig, ax = plt.subplots(figsize=(10, 7))
    y = np.arange(len(items))
    vals = [c for _, c in items]
    colors = [CLADE_COLORS.get(cl, "#64748B") for (cl, _), _ in items]
    labels = [org for (_, org), _ in items]
    ax.barh(y, vals, color=colors, alpha=0.9)
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=8)
    ax.invert_yaxis()
    ax.set_xlabel("Present / kept innexin loci")
    ax.set_title("Top species by recovered innexin loci (discovery + curator)")
    handles = [Patch(color=CLADE_COLORS[c], label=c) for c in sorted({cl for (cl, _), _ in items}) if c in CLADE_COLORS]
    ax.legend(handles=handles, fontsize=7, frameon=False, loc="lower right")
    _save(fig, out)


def plot_subfamily_by_clade(rows: list[dict], out: Path) -> None:
    sfs = ["SF1_shakB", "SF2_Inx1_ogre", "SF3_Inx3_Inx7_nematode", "SF4_Inx2_expansion", "unassigned"]
    labels = ["SF1 shakB", "SF2 Inx1/ogre", "SF3 Inx3/7+nematode", "SF4 Inx2/4–6", "unassigned"]
    colors = ["#1D4ED8", "#0F766E", "#B45309", "#0369A1", "#94A3B8"]
    clades = sorted({r["clade"] for r in rows}, key=lambda c: (-sum(1 for r in rows if r["clade"] == c), c))
    fig, ax = plt.subplots(figsize=(11, 5.5))
    x = np.arange(len(clades))
    bottom = np.zeros(len(clades))
    for sf, lab, col in zip(sfs, labels, colors):
        vals = np.array([sum(1 for r in rows if r["clade"] == c and r["subfamily"] == sf) for c in clades], dtype=float)
        ax.bar(x, vals, bottom=bottom, label=lab, color=col, width=0.72)
        bottom += vals
    ax.set_xticks(x)
    ax.set_xticklabels(clades, rotation=25, ha="right")
    ax.set_ylabel("Loci")
    ax.set_title("Tree-guided subfamily (via best-hit) across clades")
    ax.legend(frameon=False, fontsize=8)
    _save(fig, out)


def plot_inx_family_radar_like(rows: list[dict], out: Path) -> None:
    """Grouped bars: each classical Dmel type across major clades."""
    focus = ["Inx1/ogre", "Inx2", "Inx3", "Inx4/zpg", "Inx5", "Inx6", "Inx7", "shakB"]
    clades = ["Reference insect", "Insecta", "Nematoda", "Reference nematode", "Rotifera", "Mollusca", "Arthropoda (non-insect)"]
    clades = [c for c in clades if any(r["clade"] == c for r in rows)]
    fig, ax = plt.subplots(figsize=(12, 5.5))
    x = np.arange(len(focus))
    width = 0.11
    for i, clade in enumerate(clades):
        vals = [sum(1 for r in rows if r["reference_type"] == t and r["clade"] == clade) for t in focus]
        ax.bar(x + (i - len(clades) / 2) * width, vals, width=width, label=clade, color=CLADE_COLORS.get(clade, "#64748B"), alpha=0.9)
    ax.set_xticks(x)
    ax.set_xticklabels(focus, rotation=25, ha="right")
    ax.set_ylabel("Locus count")
    ax.set_title("Classical innexin labels (Inx1–7, shakB) across clades")
    ax.legend(fontsize=7, ncol=2, frameon=False)
    _save(fig, out)


def plot_identity_vs_length(rows: list[dict], out: Path) -> None:
    disc = [r for r in rows if r["source"] != "reference_db" and r["reference_identity"] > 0 and r["protein_length"] > 0]
    if len(disc) < 5:
        return
    fig, ax = plt.subplots(figsize=(8, 5.5))
    for clade in sorted({r["clade"] for r in disc}):
        pts = [r for r in disc if r["clade"] == clade]
        ax.scatter(
            [r["protein_length"] for r in pts],
            [r["reference_identity"] * 100 for r in pts],
            s=36,
            alpha=0.75,
            label=clade,
            c=CLADE_COLORS.get(clade, "#64748B"),
            edgecolors="white",
            linewidths=0.4,
        )
    ax.set_xlabel("Protein length (aa)")
    ax.set_ylabel("Best-hit identity (%)")
    ax.set_title("Discovery/curator loci: length vs identity by clade")
    ax.legend(fontsize=7, frameon=False)
    _save(fig, out)


def plot_source_donut(rows: list[dict], out: Path) -> None:
    counts = Counter(r["source"] for r in rows)
    labels = {"reference_db": "Reference DB", "discovery_old": "Old discovery", "curator_new": "Curator new"}
    colors = {"reference_db": "#115E59", "discovery_old": "#64748B", "curator_new": "#0EA5E9"}
    keys = [k for k in ("reference_db", "discovery_old", "curator_new") if counts.get(k)]
    fig, ax = plt.subplots(figsize=(6, 6))
    ax.pie(
        [counts[k] for k in keys],
        labels=[f"{labels[k]}\n({counts[k]})" for k in keys],
        colors=[colors[k] for k in keys],
        wedgeprops=dict(width=0.45),
        startangle=90,
    )
    ax.set_title("Dataset composition")
    _save(fig, out)


def write_result_tables(rows: list[dict]) -> None:
    write_csv(OUT / "all_innexin_loci_annotated.csv", rows)
    # clade summary
    clade_rows = []
    for clade in sorted({r["clade"] for r in rows}):
        sub = [r for r in rows if r["clade"] == clade]
        ids = [r["reference_identity"] * 100 for r in sub if r["reference_identity"] > 0 and r["source"] != "reference_db"]
        lens = [r["protein_length"] for r in sub if r["protein_length"] > 0]
        clade_rows.append(
            {
                "clade": clade,
                "n_loci": len(sub),
                "n_species": len({r["species_dir"] for r in sub}),
                "n_reference": sum(1 for r in sub if r["source"] == "reference_db"),
                "n_discovery_old": sum(1 for r in sub if r["source"] == "discovery_old"),
                "n_curator_new": sum(1 for r in sub if r["source"] == "curator_new"),
                "median_length_aa": f"{np.median(lens):.0f}" if lens else "",
                "median_identity_pct": f"{np.median(ids):.1f}" if ids else "",
                "top_reference_types": "; ".join(f"{k}:{v}" for k, v in Counter(r["reference_type"] for r in sub).most_common(5)),
            }
        )
    write_csv(OUT / "clade_summary.csv", clade_rows)

    type_rows = []
    for t in sorted({r["reference_type"] for r in rows}):
        sub = [r for r in rows if r["reference_type"] == t]
        type_rows.append(
            {
                "reference_type": t,
                "subfamily": SUBFAMILY_FROM_TYPE.get(t, "unassigned"),
                "n_loci": len(sub),
                "n_clades": len({r["clade"] for r in sub}),
                "clades": "; ".join(sorted({r["clade"] for r in sub})),
                "n_curator": sum(1 for r in sub if r["source"] == "curator_new"),
                "n_reference": sum(1 for r in sub if r["source"] == "reference_db"),
            }
        )
    write_csv(OUT / "reference_type_summary.csv", type_rows)

    # pairwise clade comparisons (median identity / length)
    pair_rows = []
    clades = sorted({r["clade"] for r in rows if r["source"] != "reference_db"})
    for i, a in enumerate(clades):
        for b in clades[i + 1 :]:
            la = [r["protein_length"] for r in rows if r["clade"] == a and r["protein_length"] > 0]
            lb = [r["protein_length"] for r in rows if r["clade"] == b and r["protein_length"] > 0]
            ia = [r["reference_identity"] * 100 for r in rows if r["clade"] == a and r["reference_identity"] > 0]
            ib = [r["reference_identity"] * 100 for r in rows if r["clade"] == b and r["reference_identity"] > 0]
            if not la or not lb:
                continue
            pair_rows.append(
                {
                    "clade_a": a,
                    "clade_b": b,
                    "median_length_a": f"{np.median(la):.0f}",
                    "median_length_b": f"{np.median(lb):.0f}",
                    "delta_median_length": f"{np.median(la) - np.median(lb):.0f}",
                    "median_identity_a": f"{np.median(ia):.1f}" if ia else "",
                    "median_identity_b": f"{np.median(ib):.1f}" if ib else "",
                    "n_a": len(la),
                    "n_b": len(lb),
                }
            )
    write_csv(OUT / "pairwise_clade_comparisons.csv", pair_rows)


def key_findings(rows: list[dict]) -> list[str]:
    findings = []
    disc = [r for r in rows if r["source"] != "reference_db"]
    curator = [r for r in rows if r["source"] == "curator_new"]
    findings.append(
        f"Dataset integrates {sum(1 for r in rows if r['source']=='reference_db')} reference proteins, "
        f"{sum(1 for r in rows if r['source']=='discovery_old')} old-discovery loci and "
        f"{len(curator)} curator-present loci across {len({r['species_dir'] for r in rows})} species."
    )
    # dominant type among curator
    if curator:
        top_t, n_t = Counter(r["reference_type"] for r in curator).most_common(1)[0]
        findings.append(
            f"Curator recoveries most often match {top_t} references ({n_t}/{len(curator)} loci) — "
            "names like Inx2 in nematodes are best-hit labels, not proven orthologs."
        )
    # SF3 enrichment outside insects
    non_insect = [r for r in disc if r["clade"] not in ("Insecta", "Reference insect")]
    sf3 = sum(1 for r in non_insect if r["subfamily"] == "SF3_Inx3_Inx7_nematode")
    if non_insect:
        findings.append(
            f"Outside insects, {sf3}/{len(non_insect)} kept loci map (via best-hit) to SF3 "
            "(Inx3/Inx7 + nematode-like references), consistent with the tree grouping."
        )
    # clade with most loci
    if disc:
        top_c, n_c = Counter(r["clade"] for r in disc).most_common(1)[0]
        findings.append(f"Among discovery+curator loci, {top_c} contributes the most ({n_c} loci).")
    # length contrast insects vs rotifers
    def med_len(clade: str) -> float | None:
        vals = [r["protein_length"] for r in rows if r["clade"] == clade and r["protein_length"] > 50]
        return float(np.median(vals)) if vals else None

    for a, b in (("Rotifera", "Mollusca"), ("Arthropoda (non-insect)", "Rotifera"), ("Reference insect", "Reference nematode")):
        ma, mb = med_len(a), med_len(b)
        if ma and mb:
            findings.append(f"Median protein length: {a} {ma:.0f} aa vs {b} {mb:.0f} aa.")
    # classical types present in reference
    ref_types = sorted({r["reference_type"] for r in rows if r["source"] == "reference_db"})
    findings.append("Reference DB covers classical labels: " + ", ".join(ref_types) + ".")
    findings.append(
        "Note: flat median ~1 exon in Ctenophora / non-insect Arthropoda is often a TSA / ab initio "
        "gene model (continuous ORF as one exon), not intron-free innexins."
    )
    return findings


def build_html(rows: list[dict], findings: list[str], fig_names: list[tuple[str, str]]) -> None:
    stats = summarize(rows)
    finding_html = "".join(f"<li>{html.escape(f)}</li>" for f in findings)
    cards = []
    for fname, caption in fig_names:
        if not (FIGS / fname).exists():
            continue
        cards.append(
            f"""
<article class="card">
  <img src="figures/{html.escape(fname)}" alt="{html.escape(caption)}" loading="lazy" onclick="openModal(this.src)">
  <h3>{html.escape(caption)}</h3>
</article>"""
        )

    # compact clade table
    clade_table = ["<tr><th>Clade</th><th>Loci</th><th>Species</th><th>Ref</th><th>Old</th><th>Curator</th></tr>"]
    for clade, n in stats["by_clade"].most_common():
        sub = [r for r in rows if r["clade"] == clade]
        clade_table.append(
            "<tr>"
            f"<td>{html.escape(clade)}</td><td>{n}</td>"
            f"<td>{len({r['species_dir'] for r in sub})}</td>"
            f"<td>{sum(1 for r in sub if r['source']=='reference_db')}</td>"
            f"<td>{sum(1 for r in sub if r['source']=='discovery_old')}</td>"
            f"<td>{sum(1 for r in sub if r['source']=='curator_new')}</td>"
            "</tr>"
        )

    type_table = ["<tr><th>Reference type</th><th>Subfamily</th><th>Loci</th><th>Clades</th></tr>"]
    for t, n in stats["by_type"].most_common():
        type_table.append(
            f"<tr><td>{html.escape(t)}</td><td>{html.escape(SUBFAMILY_FROM_TYPE.get(t,'unassigned'))}</td>"
            f"<td>{n}</td><td>{len({r['clade'] for r in rows if r['reference_type']==t})}</td></tr>"
        )

    page = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="utf-8">
<meta name="viewport" content="width=device-width, initial-scale=1">
<title>Innexin Clade × Type Comparison</title>
<link rel="preconnect" href="https://fonts.googleapis.com">
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600;9..144,700&family=Sora:wght@400;550;700&display=swap" rel="stylesheet">
<style>
:root {{
  --ink:#102a33; --muted:#4d6570; --sea:#0f766e; --deep:#0b3b40; --sand:#d4a017; --bg:#eef6f4;
}}
*{{box-sizing:border-box}}
body{{margin:0;font-family:Sora,sans-serif;color:var(--ink);background:
 radial-gradient(900px 420px at 0% 0%, rgba(14,165,233,.12), transparent 55%),
 radial-gradient(800px 380px at 100% 10%, rgba(15,118,110,.14), transparent 50%),
 linear-gradient(180deg,#e7f3f0,var(--bg));line-height:1.55}}
.hero{{min-height:78vh;display:grid;align-items:end;padding:clamp(1.4rem,4vw,3rem);color:#f4fffb;
 background:linear-gradient(120deg,rgba(11,59,64,.94),rgba(15,118,110,.78) 55%,rgba(180,83,9,.55)),
 radial-gradient(circle at 75% 25%,#1d8a7a,#0b3b40 60%)}}
.brand{{font-family:Fraunces,serif;font-size:clamp(2.4rem,6.5vw,4.6rem);margin:0 0 .7rem;line-height:.95;letter-spacing:-.03em}}
.hero h1{{font-size:clamp(1.05rem,2vw,1.3rem);font-weight:500;max-width:40ch;margin:0 0 .7rem}}
.hero p{{max-width:52ch;opacity:.92;margin:0 0 1.4rem}}
.cta a{{display:inline-block;text-decoration:none;background:#f4fffb;color:var(--deep);padding:.8rem 1.15rem;border-radius:999px;font-weight:700;margin-right:.5rem}}
.cta a.ghost{{background:transparent;color:#f4fffb;border:1px solid rgba(244,255,251,.45)}}
.wrap{{width:min(1160px,calc(100% - 2rem));margin:0 auto;padding:2rem 0 4rem}}
.stats{{display:grid;grid-template-columns:repeat(4,1fr);gap:.8rem;margin-top:-2rem;position:relative;z-index:2}}
.stat{{background:rgba(255,255,255,.9);border:1px solid rgba(16,42,51,.08);border-radius:16px;padding:1rem;box-shadow:0 10px 28px rgba(16,42,51,.06)}}
.stat b{{display:block;font-family:Fraunces,serif;font-size:1.8rem;color:var(--deep)}}
.stat span{{color:var(--muted);font-size:.84rem}}
section{{margin:2.4rem 0}}
h2{{font-family:Fraunces,serif;font-size:1.7rem;margin:0 0 .6rem}}
.lead{{color:var(--muted);max-width:70ch}}
.findings{{background:#fff;border-radius:16px;padding:1.1rem 1.3rem;border-left:4px solid var(--sand);box-shadow:0 8px 24px rgba(16,42,51,.05)}}
.findings ol{{margin:.4rem 0 0;padding-left:1.2rem}}
.findings li{{margin:.35rem 0;color:var(--muted)}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(320px,1fr));gap:1rem}}
.card{{background:rgba(255,255,255,.92);border:1px solid rgba(16,42,51,.08);border-radius:16px;padding:.85rem;box-shadow:0 8px 22px rgba(16,42,51,.05)}}
.card img{{width:100%;border-radius:10px;cursor:zoom-in;background:#f8fafc}}
.card h3{{font-size:.95rem;margin:.7rem 0 0;font-weight:600}}
table{{width:100%;border-collapse:collapse;background:#fff;border-radius:12px;overflow:hidden;box-shadow:0 6px 18px rgba(16,42,51,.05);font-size:.88rem}}
th,td{{padding:.55rem .7rem;border-bottom:1px solid #e2e8f0;text-align:left}}
th{{background:#0b3b40;color:#f4fffb;font-weight:600}}
.links a{{margin-right:1rem}}
.modal{{display:none;position:fixed;inset:0;background:rgba(0,0,0,.88);z-index:20;align-items:center;justify-content:center;padding:1.5rem;cursor:zoom-out}}
.modal.open{{display:flex}}
.modal img{{max-width:95vw;max-height:90vh;border-radius:8px}}
@media(max-width:800px){{.stats{{grid-template-columns:1fr 1fr}}.hero{{min-height:68vh}}}}
</style>
</head>
<body>
<header class="hero">
  <div>
    <p class="brand">Clade × Innexin</p>
    <h1>Compare every clade and every classical innexin label in one place</h1>
    <p>Reference database proteins, older discovery candidates and curator-recovered loci — sliced by clade, best-hit type (Inx1–7, shakB, unc-7/9…) and tree-guided subfamily.</p>
    <div class="cta">
      <a href="#figures">See diagrams</a>
      <a class="ghost" href="../new_species_gallery/index.html">New species gallery</a>
    </div>
  </div>
</header>
<main class="wrap">
  <div class="stats">
    <div class="stat"><b>{stats['n_loci']}</b><span>annotated innexin loci</span></div>
    <div class="stat"><b>{stats['n_species']}</b><span>species represented</span></div>
    <div class="stat"><b>{len(stats['by_clade'])}</b><span>clade bins</span></div>
    <div class="stat"><b>{len(stats['by_type'])}</b><span>reference-type labels</span></div>
  </div>

  <section>
    <h2>Key results</h2>
    <div class="findings"><ol>{finding_html}</ol></div>
    <p class="lead" style="margin-top:1rem">
      “Inx2-like” for a mollusc/rotifer means the locus’s <em>best reference hit</em> was an Inx2/unc-9/… protein.
      Orthology still needs the phylogeny (SF1–SF4), not the gene name alone.
    </p>
  </section>

  <section>
    <h2>Clade overview</h2>
    <table>{''.join(clade_table)}</table>
  </section>

  <section>
    <h2>Innexin labels (Inx1–7, shakB, nematode genes)</h2>
    <table>{''.join(type_table)}</table>
  </section>

  <section id="figures">
    <h2>Diagrams</h2>
    <p class="lead">Every comparison plot generated from the combined reference + discovery + curator table.</p>
    <div class="grid">
      {''.join(cards)}
    </div>
  </section>

  <section class="links">
    <h2>Data tables</h2>
    <p>
      <a href="all_innexin_loci_annotated.csv">all_innexin_loci_annotated.csv</a>
      <a href="clade_summary.csv">clade_summary.csv</a>
      <a href="reference_type_summary.csv">reference_type_summary.csv</a>
      <a href="pairwise_clade_comparisons.csv">pairwise_clade_comparisons.csv</a>
    </p>
    <p>
      <a href="../phylogenetic_story/index.html">Phylogeny story</a>
      <a href="../innexin_subfamilies/index.html">Subfamilies</a>
      <a href="../family_comparison/index.html">vs connexin</a>
      <a href="../showcase/index.html">Full showcase</a>
    </p>
  </section>
</main>
<div id="modal" class="modal" onclick="this.classList.remove('open')"><img id="modal-img" alt=""></div>
<script>
function openModal(src){{document.getElementById('modal-img').src=src;document.getElementById('modal').classList.add('open')}}
document.addEventListener('keydown',e=>{{if(e.key==='Escape')document.getElementById('modal').classList.remove('open')}});
</script>
</body>
</html>
"""
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "index.html").write_text(page, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)

    rows = load_reference_rows() + load_candidate_rows()
    if not rows:
        raise SystemExit("No innexin rows found")

    write_result_tables(rows)

    fig_plan = [
        ("01_loci_by_clade_source.png", "Loci by clade, stacked by source"),
        ("02_heatmap_clade_x_type.png", "Heatmap: clade × reference type"),
        ("03_type_mix_per_clade.png", "Reference-type mix within clades"),
        ("04_classical_inx_across_clades.png", "Inx1–7 & shakB across clades"),
        ("05_subfamily_by_clade.png", "Tree subfamilies SF1–SF4 by clade"),
        ("06_length_by_clade.png", "Protein length by clade"),
        ("07_length_by_type.png", "Protein length by innexin label"),
        ("08_identity_by_clade.png", "Best-hit identity (%) by clade"),
        ("09_identity_by_type.png", "Best-hit identity (%) by innexin label"),
        ("10_species_richness.png", "Top species by locus count"),
        ("11_identity_vs_length.png", "Length vs identity scatter"),
        ("12_dataset_composition.png", "Dataset composition"),
        ("13_exon_by_clade.png", "Exon count by clade"),
        ("14_tm_by_type.png", "TM helices by innexin label (where annotated)"),
    ]

    plot_loci_by_clade(rows, FIGS / "01_loci_by_clade_source.png")
    plot_heatmap_clade_type(rows, FIGS / "02_heatmap_clade_x_type.png")
    plot_type_composition_per_clade(rows, FIGS / "03_type_mix_per_clade.png")
    plot_inx_family_radar_like(rows, FIGS / "04_classical_inx_across_clades.png")
    plot_subfamily_by_clade(rows, FIGS / "05_subfamily_by_clade.png")
    _boxplot_by(rows, "clade", "protein_length", FIGS / "06_length_by_clade.png", "Protein length by clade", "Length (aa)", min_n=2)
    _boxplot_by(rows, "reference_type", "protein_length", FIGS / "07_length_by_type.png", "Protein length by innexin label", "Length (aa)", min_n=2)
    _boxplot_by(
        [r for r in rows if r["source"] != "reference_db"],
        "clade",
        "reference_identity",
        FIGS / "08_identity_by_clade.png",
        "Best-hit identity by clade (discovery + curator)",
        "Identity (%)",
        min_n=2,
    )
    _boxplot_by(
        [r for r in rows if r["source"] != "reference_db"],
        "reference_type",
        "reference_identity",
        FIGS / "09_identity_by_type.png",
        "Best-hit identity by innexin label",
        "Identity (%)",
        min_n=2,
    )
    plot_species_richness(rows, FIGS / "10_species_richness.png")
    plot_identity_vs_length(rows, FIGS / "11_identity_vs_length.png")
    plot_source_donut(rows, FIGS / "12_dataset_composition.png")
    _boxplot_by(
        [r for r in rows if r["exon_count"]],
        "clade",
        "exon_count",
        FIGS / "13_exon_by_clade.png",
        "Exon count by clade (1-exon medians often = TSA/ab initio models)",
        "Exons",
        min_n=2,
    )
    _boxplot_by(
        [r for r in rows if r["tm_helix_count"]],
        "reference_type",
        "tm_helix_count",
        FIGS / "14_tm_by_type.png",
        "Predicted TM helices by innexin label",
        "TM helices",
        min_n=2,
    )

    findings = key_findings(rows)
    (OUT / "key_findings.md").write_text("# Key findings\n\n" + "\n".join(f"- {f}" for f in findings) + "\n", encoding="utf-8")
    build_html(rows, findings, fig_plan)
    print(f"Wrote {OUT}/index.html")
    print(f"Loci={len(rows)} species={len({r['species_dir'] for r in rows})} figures={len(list(FIGS.glob('*.png')))}")


if __name__ == "__main__":
    main()
