#!/usr/bin/env python3
"""Protein-level (MMseqs2) + embedding-level similarity across innexin groups.

Approaches:
  1. MMseqs2 all-vs-all search (sequence identity / bit score)
  2. MMseqs2 cluster threshold sweep (0.20–0.70) — clade-thinking style
  3. AA-composition embeddings + cosine similarity
  4. 3-mer frequency embeddings (PCA/SVD) + cosine similarity

Outputs → project/results/innexin_similarity/
"""

from __future__ import annotations

import csv
import html
import math
import re
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from Bio.Seq import Seq

from pipeline.batch_innexin_curator_probe import parse_miniprot_models, translate_model
from pipeline.common import METADATA_DIR, PROJECT_ROOT, RESULTS_DIR, write_csv
from pipeline.innexin_clade_comparison_builder import (
    SUBFAMILY_FROM_TYPE,
    classify_reference_type,
    infer_clade,
)

OUT = RESULTS_DIR / "innexin_similarity"
FIGS = OUT / "figures"
WORK = OUT / "work"
FASTA = OUT / "all_innexins.fasta"
META = OUT / "sequence_metadata.csv"

MMSEQS = shutil.which("mmseqs") or str(
    Path.home() / "miniconda3/envs/synvoy_env/bin/mmseqs"
)
SAMTOOLS = shutil.which("samtools") or str(
    Path.home() / "miniconda3/envs/synvoy_env/bin/samtools"
)

REF_ROOT = PROJECT_ROOT / "project/data/references/innexins"
CURATOR = RESULTS_DIR / "gene_curator_probe"
DISCOVERY = RESULTS_DIR / "innexin_discovery"
PHYLO_FASTA = RESULTS_DIR / "phylogeny" / "innexin_reference_panel.fasta"
SUMMARY = CURATOR / "innexin_search_summary.csv"

# Tip / gene → tree subfamily (reference anchors)
GENE_TO_SF = {
    "ogre": "SF2_Inx1_ogre",
    "inx1": "SF2_Inx1_ogre",
    "inx-1": "SF2_Inx1_ogre",
    "inx2": "SF4_Inx2_expansion",
    "inx-2": "SF3_Inx3_Inx7_nematode",  # Cele inx-2 is SF3; insect Inx2 handled separately
    "inx3": "SF3_Inx3_Inx7_nematode",
    "inx-3": "SF3_Inx3_Inx7_nematode",
    "zpg": "SF4_Inx2_expansion",
    "inx4": "SF4_Inx2_expansion",
    "inx5": "SF4_Inx2_expansion",
    "inx6": "SF4_Inx2_expansion",
    "inx7": "SF3_Inx3_Inx7_nematode",
    "inx-7": "SF3_Inx3_Inx7_nematode",
    "shakb": "SF1_shakB",
    "unc-7": "SF3_Inx3_Inx7_nematode",
    "unc-9": "SF3_Inx3_Inx7_nematode",
    "eat-5": "SF3_Inx3_Inx7_nematode",
}

SF_COLORS = {
    "SF1_shakB": "#1D4ED8",
    "SF2_Inx1_ogre": "#0F766E",
    "SF3_Inx3_Inx7_nematode": "#B45309",
    "SF4_Inx2_expansion": "#0369A1",
    "unassigned": "#94A3B8",
}

THRESHOLDS = [0.20, 0.30, 0.40, 0.50, 0.60, 0.70]

plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _save(fig: plt.Figure, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def _clean_aa(seq: str) -> str:
    s = re.sub(r"[^A-Za-z*]", "", seq).upper().split("*")[0]
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", s)


def _sf_from_gene(gene: str, organism: str = "") -> str:
    g = (gene or "").lower().strip()
    org = (organism or "").lower()
    if g in ("inx2", "inx-2") and ("drosophila" in org or "schistocerca" in org or "anopheles" in org):
        return "SF4_Inx2_expansion"
    if g in ("inx2", "inx-2") and "caenorhabditis" in org:
        return "SF3_Inx3_Inx7_nematode"
    if g.startswith("inx-") and "caenorhabditis" in org:
        return "SF3_Inx3_Inx7_nematode"
    return GENE_TO_SF.get(g, SUBFAMILY_FROM_TYPE.get(classify_reference_type(gene, gene), "unassigned"))


def _write_fa(handle, seq_id: str, seq: str) -> None:
    seq = _clean_aa(seq)
    if len(seq) < 80:
        return
    handle.write(f">{seq_id}\n")
    for i in range(0, len(seq), 80):
        handle.write(seq[i : i + 80] + "\n")


def _read_fasta(path: Path) -> list[tuple[str, str]]:
    out: list[tuple[str, str]] = []
    if not path.exists():
        return out
    hdr, chunks = None, []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            if hdr is not None:
                out.append((hdr, "".join(chunks)))
            hdr = line[1:].strip()
            chunks = []
        else:
            chunks.append(line.strip())
    if hdr is not None:
        out.append((hdr, "".join(chunks)))
    return out


def collect_sequences() -> list[dict]:
    """Build annotated protein set from reference, phylogeny, discovery, curator."""
    rows: list[dict] = []
    seen_seq: set[str] = set()

    def add(seq_id: str, seq: str, meta: dict) -> None:
        aa = _clean_aa(seq)
        if len(aa) < 80:
            return
        key = f"{meta.get('species_dir')}|{aa[:40]}|{len(aa)}"
        if key in seen_seq:
            return
        seen_seq.add(key)
        meta = dict(meta)
        meta["seq_id"] = seq_id
        meta["length"] = len(aa)
        meta["seq"] = aa
        rows.append(meta)

    # 1) Reference FASTAs
    if REF_ROOT.exists():
        for fasta in sorted(REF_ROOT.rglob("*.fasta")):
            organism = fasta.parent.name.replace("_", " ")
            species_dir = fasta.parent.name
            gene = fasta.stem.split("__")[0]
            acc = fasta.stem.split("__")[1] if "__" in fasta.stem else ""
            for _, seq in _read_fasta(fasta):
                sid = f"ref|{species_dir}|{gene}|{acc or 'na'}"
                add(
                    sid,
                    seq,
                    {
                        "source": "reference_db",
                        "species_dir": species_dir,
                        "organism": organism,
                        "clade": infer_clade(species_dir, organism),
                        "gene_label": gene,
                        "reference_type": classify_reference_type(fasta.stem, gene),
                        "subfamily": _sf_from_gene(gene, organism),
                    },
                )

    # 2) Phylogeny panel (may overlap refs; dedupe by sequence prefix)
    for hdr, seq in _read_fasta(PHYLO_FASTA):
        # headers like sp|Q9V427|INX2_DROME_Innexin or Species|candidate...
        gene = hdr
        organism = ""
        species_dir = "phylo_panel"
        if "|" in hdr:
            parts = hdr.split("|")
            if parts[0] == "sp" and len(parts) >= 3:
                gene = parts[2].split("_")[0]
                if "DROME" in parts[2]:
                    organism = "Drosophila melanogaster"
                    species_dir = "Drosophila_melanogaster"
                elif "CAEEL" in parts[2]:
                    organism = "Caenorhabditis elegans"
                    species_dir = "Caenorhabditis_elegans"
                elif "CAEBR" in parts[2]:
                    organism = "Caenorhabditis briggsae"
                    species_dir = "Caenorhabditis_briggsae"
                elif "ANOGA" in parts[2]:
                    organism = "Anopheles gambiae"
                    species_dir = "Anopheles_gambiae"
                elif "AEDAE" in parts[2]:
                    organism = "Aedes aegypti"
                    species_dir = "Aedes_aegypti"
                elif "SCHAM" in parts[2]:
                    organism = "Schistocerca americana"
                    species_dir = "Schistocerca_americana"
            else:
                species_dir = parts[0].replace(" ", "_")
                organism = parts[0].replace("_", " ")
                gene = parts[1] if len(parts) > 1 else hdr
        sid = f"phylo|{re.sub(r'[^A-Za-z0-9_.-]+', '_', hdr)[:80]}"
        add(
            sid,
            seq,
            {
                "source": "phylo_panel",
                "species_dir": species_dir,
                "organism": organism or species_dir.replace("_", " "),
                "clade": infer_clade(species_dir, organism),
                "gene_label": gene,
                "reference_type": classify_reference_type(hdr, gene),
                "subfamily": _sf_from_gene(gene, organism),
            },
        )

    # 3) Discovery candidate FASTAs
    if DISCOVERY.exists():
        for fasta in sorted(DISCOVERY.glob("*/candidate_innexins.fasta")):
            species_dir = fasta.parent.name
            organism = species_dir.replace("_", " ")
            for i, (hdr, seq) in enumerate(_read_fasta(fasta), 1):
                sid = f"disc|{species_dir}|{re.sub(r'[^A-Za-z0-9_.-]+', '_', hdr)[:60] or f'c{i}'}"
                rtype = classify_reference_type(hdr, hdr)
                add(
                    sid,
                    seq,
                    {
                        "source": "discovery_old",
                        "species_dir": species_dir,
                        "organism": organism,
                        "clade": infer_clade(species_dir, organism),
                        "gene_label": hdr.split()[0][:40],
                        "reference_type": rtype,
                        "subfamily": SUBFAMILY_FROM_TYPE.get(rtype, "unassigned"),
                    },
                )

    # 4) Curator present products (retranslate from miniprot GFF; present loci only)
    if SUMMARY.exists():
        with SUMMARY.open(newline="", encoding="utf-8") as fh:
            for rec in csv.DictReader(fh):
                if int(rec.get("present_loci") or 0) <= 0:
                    continue
                species = rec["species"]
                genome = Path(rec["genome"])
                if not genome.is_absolute():
                    genome = PROJECT_ROOT / genome
                gff = CURATOR / species / "miniprot.gff"
                hits_csv = CURATOR / species / "hits.csv"
                if not genome.exists() or not gff.exists() or not hits_csv.exists():
                    continue
                present_windows: set[tuple[str, int, int]] = set()
                with hits_csv.open(newline="", encoding="utf-8") as hf:
                    for h in csv.DictReader(hf):
                        if (h.get("verdict") or "") != "present":
                            continue
                        present_windows.add((h["seqid"], int(h["start"]), int(h["end"])))
                models = parse_miniprot_models(gff)
                n = 0
                for model in models:
                    key = (str(model.get("seqid")), int(model.get("start") or 0), int(model.get("end") or 0))
                    if key not in present_windows:
                        continue
                    try:
                        prot = translate_model(genome, model)
                    except Exception:
                        continue
                    aa = _clean_aa(prot)
                    if len(aa) < 180:
                        continue
                    identity = float(model.get("identity") or 0)
                    n += 1
                    target = str(model.get("target") or "")
                    rtype = classify_reference_type(target, target)
                    sid = (
                        f"cur|{species}|{model.get('seqid')}_"
                        f"{model.get('start')}_{model.get('end')}|p{n}"
                    )
                    add(
                        sid,
                        aa,
                        {
                            "source": "curator_new",
                            "species_dir": species,
                            "organism": species.replace("_", " "),
                            "clade": infer_clade(species, species.replace("_", " ")),
                            "gene_label": target.split()[0] if target else f"curator_{n}",
                            "reference_type": rtype,
                            "subfamily": SUBFAMILY_FROM_TYPE.get(rtype, "unassigned"),
                            "identity_to_ref": identity,
                        },
                    )

    return rows


def write_fasta_and_meta(rows: list[dict]) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    with FASTA.open("w", encoding="utf-8") as fh:
        for r in rows:
            _write_fa(fh, r["seq_id"], r["seq"])
    fields = [
        "seq_id",
        "source",
        "species_dir",
        "organism",
        "clade",
        "gene_label",
        "reference_type",
        "subfamily",
        "length",
        "identity_to_ref",
    ]
    meta_rows = []
    for r in rows:
        meta_rows.append({k: r.get(k, "") for k in fields})
    write_csv(META, meta_rows, fieldnames=fields)


def run_mmseqs(rows: list[dict]) -> tuple[Path, dict[float, Path]]:
    WORK.mkdir(parents=True, exist_ok=True)
    db = WORK / "seqdb"
    tmp = WORK / "tmp"
    tmp.mkdir(exist_ok=True)
    if db.with_suffix(".dbtype").exists() or (WORK / "seqdb").exists():
        # recreate cleanly
        for p in WORK.glob("seqdb*"):
            if p.is_dir():
                shutil.rmtree(p)
            else:
                p.unlink()
    subprocess.run([MMSEQS, "createdb", str(FASTA), str(db)], check=True)

    # all-vs-all search
    search_db = WORK / "allvsall"
    search_tsv = OUT / "mmseqs_allvsall.m8"
    if search_db.with_suffix(".dbtype").exists() or (WORK / "allvsall.index").exists():
        for p in WORK.glob("allvsall*"):
            if p.is_file():
                p.unlink()
    subprocess.run(
        [
            MMSEQS,
            "search",
            str(db),
            str(db),
            str(search_db),
            str(tmp),
            "-s",
            "7.5",
            "--max-seqs",
            "1000",
            "-e",
            "1e-3",
            "-v",
            "1",
        ],
        check=True,
    )
    subprocess.run(
        [
            MMSEQS,
            "convertalis",
            str(db),
            str(db),
            str(search_db),
            str(search_tsv),
            "--format-output",
            "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
            "-v",
            "1",
        ],
        check=True,
    )

    cluster_tsvs: dict[float, Path] = {}
    for thr in THRESHOLDS:
        tag = f"clust_{int(thr * 100):02d}"
        clust = WORK / tag
        tsv = OUT / f"mmseqs_clusters_id{int(thr * 100):02d}.tsv"
        for p in WORK.glob(f"{tag}*"):
            if p.is_file():
                p.unlink()
        subprocess.run(
            [
                MMSEQS,
                "cluster",
                str(db),
                str(clust),
                str(tmp),
                "--min-seq-id",
                str(thr),
                "-c",
                "0.8",
                "--cov-mode",
                "0",
                "-v",
                "1",
            ],
            check=True,
        )
        subprocess.run(
            [MMSEQS, "createtsv", str(db), str(db), str(clust), str(tsv), "-v", "1"],
            check=True,
        )
        cluster_tsvs[thr] = tsv

    return search_tsv, cluster_tsvs


def load_meta_map(rows: list[dict]) -> dict[str, dict]:
    return {r["seq_id"]: r for r in rows}


def summarize_clusters(cluster_tsvs: dict[float, Path], meta: dict[str, dict]) -> list[dict]:
    out = []
    for thr, path in sorted(cluster_tsvs.items()):
        groups: dict[str, list[str]] = defaultdict(list)
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            rep, mem = line.split("\t")[:2]
            groups[rep].append(mem)
        n_clust = len(groups)
        sizes = [len(v) for v in groups.values()]
        # clean-orthology style: among clusters with >=2 members, fraction species with exactly 1 member
        clean_scores = []
        sf_purity = []
        for members in groups.values():
            if len(members) < 2:
                continue
            spp = [meta.get(m, {}).get("species_dir", "?") for m in members]
            counts = Counter(spp)
            if not counts:
                continue
            clean_scores.append(sum(1 for c in counts.values() if c == 1) / len(counts))
            sfs = [meta.get(m, {}).get("subfamily", "unassigned") for m in members]
            top = Counter(sfs).most_common(1)[0][1]
            sf_purity.append(top / len(sfs))
        out.append(
            {
                "min_seq_id": thr,
                "n_clusters": n_clust,
                "n_singletons": sum(1 for s in sizes if s == 1),
                "median_cluster_size": float(np.median(sizes)) if sizes else 0,
                "max_cluster_size": max(sizes) if sizes else 0,
                "mean_clean_orthology": float(np.mean(clean_scores)) if clean_scores else float("nan"),
                "mean_subfamily_purity": float(np.mean(sf_purity)) if sf_purity else float("nan"),
                "n_multi_clusters": len(clean_scores),
            }
        )
    write_csv(OUT / "mmseqs_threshold_sweep_summary.csv", out)
    return out


def pairwise_from_m8(m8: Path, meta: dict[str, dict]) -> dict:
    """Aggregate mean identity within/between subfamilies."""
    within = defaultdict(list)
    between = defaultdict(list)
    all_pairs = []
    with m8.open(encoding="utf-8") as fh:
        for line in fh:
            q, t, pident, *_ = line.split("\t")
            if q == t:
                continue
            pid = float(pident)
            sq = meta.get(q, {}).get("subfamily", "unassigned")
            st = meta.get(t, {}).get("subfamily", "unassigned")
            all_pairs.append((q, t, pid, sq, st))
            if sq == st:
                within[sq].append(pid)
            else:
                between[frozenset((sq, st))].append(pid)
                between[("ANY_BETWEEN",)].append(pid)
    return {"within": within, "between": between, "pairs": all_pairs}


def aa_composition_matrix(rows: list[dict]) -> tuple[np.ndarray, list[str]]:
    aas = list("ACDEFGHIKLMNPQRSTVWY")
    ids = [r["seq_id"] for r in rows]
    X = np.zeros((len(rows), len(aas)), dtype=float)
    for i, r in enumerate(rows):
        seq = r["seq"]
        c = Counter(seq)
        n = max(len(seq), 1)
        for j, a in enumerate(aas):
            X[i, j] = c.get(a, 0) / n
    # L2 normalize
    norms = np.linalg.norm(X, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return X / norms, ids


def kmer_embedding_matrix(rows: list[dict], k: int = 3, n_comp: int = 20) -> tuple[np.ndarray, list[str]]:
    """Sparse-ish 3-mer counts → truncated SVD via numpy."""
    ids = [r["seq_id"] for r in rows]
    # vocabulary of observed kmers (cap)
    vocab_count: Counter[str] = Counter()
    seqs = [r["seq"] for r in rows]
    for seq in seqs:
        for i in range(max(0, len(seq) - k + 1)):
            km = seq[i : i + k]
            if "X" in km:
                continue
            vocab_count[km] += 1
    vocab = [km for km, _ in vocab_count.most_common(800)]
    idx = {km: j for j, km in enumerate(vocab)}
    X = np.zeros((len(seqs), len(vocab)), dtype=float)
    for i, seq in enumerate(seqs):
        for j in range(max(0, len(seq) - k + 1)):
            km = seq[j : j + k]
            if km in idx:
                X[i, idx[km]] += 1
        s = X[i].sum()
        if s:
            X[i] /= s
    # center + SVD
    Xc = X - X.mean(axis=0, keepdims=True)
    # economy SVD on smaller side
    try:
        u, s, vt = np.linalg.svd(Xc, full_matrices=False)
        n_comp = min(n_comp, u.shape[1])
        emb = u[:, :n_comp] * s[:n_comp]
    except np.linalg.LinAlgError:
        emb = Xc[:, :n_comp] if Xc.shape[1] >= n_comp else Xc
    norms = np.linalg.norm(emb, axis=1, keepdims=True)
    norms[norms == 0] = 1
    return emb / norms, ids


def cosine_sim(A: np.ndarray) -> np.ndarray:
    return np.clip(A @ A.T, -1.0, 1.0)


def group_similarity_table(sim: np.ndarray, ids: list[str], meta: dict[str, dict], label: str) -> list[dict]:
    sfs = sorted({meta[i]["subfamily"] for i in ids if i in meta})
    id_to_i = {s: i for i, s in enumerate(ids)}
    rows = []
    for a in sfs:
        members_a = [s for s in ids if meta.get(s, {}).get("subfamily") == a]
        for b in sfs:
            members_b = [s for s in ids if meta.get(s, {}).get("subfamily") == b]
            vals = []
            for qa in members_a:
                ia = id_to_i[qa]
                for qb in members_b:
                    if qa == qb:
                        continue
                    ib = id_to_i[qb]
                    vals.append(float(sim[ia, ib]))
            if not vals:
                continue
            rows.append(
                {
                    "method": label,
                    "group_a": a,
                    "group_b": b,
                    "n_pairs": len(vals),
                    "mean_similarity": f"{np.mean(vals):.4f}",
                    "median_similarity": f"{np.median(vals):.4f}",
                    "same_group": str(a == b).lower(),
                }
            )
    return rows


def plot_threshold_sweep(summary: list[dict]) -> None:
    if not summary:
        return
    thr = [r["min_seq_id"] for r in summary]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(thr, [r["n_clusters"] for r in summary], "o-", label="# clusters", color="#0F766E")
    ax2 = ax.twinx()
    purity = [r["mean_subfamily_purity"] for r in summary]
    clean = [r["mean_clean_orthology"] for r in summary]
    ax2.plot(thr, purity, "s--", label="subfamily purity", color="#B45309")
    ax2.plot(thr, clean, "^--", label="clean-orthology score", color="#1D4ED8")
    ax.set_xlabel("MMseqs2 min-seq-id threshold")
    ax.set_ylabel("Clusters")
    ax2.set_ylabel("Score (0–1)")
    ax.set_title("Identity-threshold sweep (clade-thinking style)")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8)
    _save(fig, FIGS / "01_mmseqs_threshold_sweep.png")


def plot_within_between_identity(pw: dict) -> None:
    within = pw["within"]
    labels = sorted(within.keys())
    data = [within[k] for k in labels if within[k]]
    labs = [k for k in labels if within[k]]
    if not data:
        return
    fig, ax = plt.subplots(figsize=(9, 5))
    bp = ax.boxplot(data, tick_labels=labs, patch_artist=True, showfliers=False)
    for patch, lab in zip(bp["boxes"], labs):
        patch.set_facecolor(SF_COLORS.get(lab, "#94A3B8"))
        patch.set_alpha(0.8)
    bet = pw["between"].get(("ANY_BETWEEN",), [])
    if bet:
        ax.axhline(np.median(bet), color="#64748B", ls="--", lw=1.2, label=f"median between-SF ({np.median(bet):.1f}%)")
        ax.legend(frameon=False)
    ax.set_ylabel("MMseqs2 % identity")
    ax.set_title("Within-subfamily sequence identity")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=25, ha="right")
    _save(fig, FIGS / "02_within_sf_identity.png")


def plot_sf_heatmap(table_rows: list[dict], method: str, out_name: str, title: str, pct: bool = False) -> None:
    sub = [r for r in table_rows if r["method"] == method]
    if not sub:
        return
    groups = sorted({r["group_a"] for r in sub} | {r["group_b"] for r in sub})
    mat = np.full((len(groups), len(groups)), np.nan)
    idx = {g: i for i, g in enumerate(groups)}
    for r in sub:
        mat[idx[r["group_a"]], idx[r["group_b"]]] = float(r["mean_similarity"])
    fig, ax = plt.subplots(figsize=(7.2, 6))
    im = ax.imshow(mat, cmap="YlGnBu", vmin=np.nanmin(mat), vmax=np.nanmax(mat))
    ax.set_xticks(range(len(groups)))
    ax.set_yticks(range(len(groups)))
    short = [g.replace("SF3_Inx3_Inx7_nematode", "SF3").replace("SF4_Inx2_expansion", "SF4").replace("SF2_Inx1_ogre", "SF2").replace("SF1_shakB", "SF1") for g in groups]
    ax.set_xticklabels(short, rotation=40, ha="right", fontsize=8)
    ax.set_yticklabels(short, fontsize=8)
    for i in range(len(groups)):
        for j in range(len(groups)):
            v = mat[i, j]
            if not math.isnan(v):
                txt = f"{v:.0f}" if pct else f"{v:.2f}"
                ax.text(j, i, txt, ha="center", va="center", fontsize=7)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    _save(fig, FIGS / out_name)


def plot_embedding_pca(emb: np.ndarray, ids: list[str], meta: dict[str, dict], title: str, out_name: str) -> None:
    if emb.shape[1] < 2:
        return
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    for sf, color in SF_COLORS.items():
        pts = [i for i, s in enumerate(ids) if meta.get(s, {}).get("subfamily") == sf]
        if not pts:
            continue
        ax.scatter(emb[pts, 0], emb[pts, 1], s=28, alpha=0.75, c=color, label=sf, edgecolors="white", linewidths=0.3)
    ax.set_xlabel("Component 1")
    ax.set_ylabel("Component 2")
    ax.set_title(title)
    ax.legend(fontsize=7, frameon=False)
    _save(fig, FIGS / out_name)


def plot_method_agreement(mm_within: dict, emb_cos: dict) -> None:
    """Scatter: mean within-SF mmseqs identity vs embedding cosine."""
    sfs = sorted(set(mm_within) & set(emb_cos))
    if len(sfs) < 2:
        return
    x = [np.mean(mm_within[s]) for s in sfs]
    y = [np.mean(emb_cos[s]) for s in sfs]
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    for sf, xv, yv in zip(sfs, x, y):
        ax.scatter([xv], [yv], s=90, c=SF_COLORS.get(sf, "#64748B"), label=sf)
        ax.annotate(sf.split("_")[0], (xv, yv), fontsize=8, xytext=(4, 4), textcoords="offset points")
    ax.set_xlabel("Mean within-SF MMseqs2 % identity")
    ax.set_ylabel("Mean within-SF embedding cosine")
    ax.set_title("Sequence vs embedding cohesion by subfamily")
    ax.legend(fontsize=7, frameon=False)
    _save(fig, FIGS / "07_seq_vs_embedding_cohesion.png")


def key_findings(sweep: list[dict], pw: dict, group_tables: list[dict]) -> list[str]:
    findings = []
    if sweep:
        best = max(sweep, key=lambda r: (r["mean_subfamily_purity"] if not math.isnan(r["mean_subfamily_purity"]) else -1, -r["n_clusters"]))
        findings.append(
            f"MMseqs2 sweep: highest mean subfamily purity {best['mean_subfamily_purity']:.2f} "
            f"at min-seq-id={best['min_seq_id']:.2f} ({best['n_clusters']} clusters)."
        )
    within = pw["within"]
    bet = pw["between"].get(("ANY_BETWEEN",), [])
    if within and bet:
        wvals = [np.median(v) for v in within.values() if v]
        findings.append(
            f"Within-subfamily median identity ~{np.mean(wvals):.1f}% vs between-subfamily "
            f"median {np.median(bet):.1f}% — groups are sequence-separable but still one family."
        )
    # embedding same-group vs different
    emb_rows = [r for r in group_tables if r["method"] == "kmer3_svd_cosine"]
    same = [float(r["mean_similarity"]) for r in emb_rows if r["same_group"] == "true"]
    diff = [float(r["mean_similarity"]) for r in emb_rows if r["same_group"] == "false"]
    if same and diff:
        findings.append(
            f"3-mer embedding cosine: same-subfamily mean {np.mean(same):.3f} vs "
            f"different-subfamily {np.mean(diff):.3f}."
        )
    aa_rows = [r for r in group_tables if r["method"] == "aa_comp_cosine"]
    same_a = [float(r["mean_similarity"]) for r in aa_rows if r["same_group"] == "true"]
    diff_a = [float(r["mean_similarity"]) for r in aa_rows if r["same_group"] == "false"]
    if same_a and diff_a:
        findings.append(
            f"AA-composition cosine is weaker at separating groups "
            f"(same {np.mean(same_a):.3f} vs diff {np.mean(diff_a):.3f}) — "
            "composition alone is not enough; k-mers / MMseqs carry more signal."
        )
    # SF3 vs SF4 contrast if present
    for method in ("mmseqs_pident", "kmer3_svd_cosine"):
        sub = [r for r in group_tables if r["method"] == method]
        cross = [
            r
            for r in sub
            if {r["group_a"], r["group_b"]} == {"SF3_Inx3_Inx7_nematode", "SF4_Inx2_expansion"}
        ]
        if cross:
            findings.append(
                f"{method}: SF3↔SF4 mean similarity {cross[0]['mean_similarity']} "
                f"(n={cross[0]['n_pairs']}) — tests insect expansion vs nematode-like clade."
            )
            break
    findings.append(
        "Note: curator subfamily labels come from best-hit mapping; embedding/MMseqs "
        "separation supports the tree groups, not final orthology."
    )
    return findings


def build_html(n_seq: int, findings: list[str], fig_plan: list[tuple[str, str]]) -> None:
    cards = []
    for fn, cap in fig_plan:
        if (FIGS / fn).exists():
            cards.append(
                f"<article class='card'><img src='figures/{fn}' alt='{html.escape(cap)}' "
                f"loading='lazy' onclick=\"openModal(this.src)\"><h3>{html.escape(cap)}</h3></article>"
            )
    finding_li = "".join(f"<li>{html.escape(f)}</li>" for f in findings)
    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Innexin group similarity — MMseqs2 & embeddings</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root{{--ink:#102a33;--muted:#4d6570;--deep:#0b3b40;--sea:#0f766e;--bg:#eef6f4}}
body{{margin:0;font-family:Sora,sans-serif;color:var(--ink);background:linear-gradient(180deg,#e7f3f0,var(--bg));line-height:1.55}}
.hero{{min-height:70vh;display:grid;align-items:end;padding:clamp(1.3rem,4vw,3rem);color:#f4fffb;
background:linear-gradient(115deg,rgba(11,59,64,.94),rgba(15,118,110,.75) 55%,rgba(3,105,161,.55))}}
.brand{{font-family:Fraunces,serif;font-size:clamp(2.2rem,6vw,4.2rem);margin:0 0 .6rem;line-height:.95}}
.hero p{{max-width:52ch;opacity:.92}}
.cta a{{display:inline-block;margin:.4rem .4rem 0 0;padding:.75rem 1.1rem;border-radius:999px;background:#f4fffb;color:var(--deep);text-decoration:none;font-weight:700}}
.wrap{{width:min(1140px,calc(100% - 2rem));margin:0 auto;padding:2rem 0 4rem}}
.stats{{display:grid;grid-template-columns:repeat(3,1fr);gap:.8rem;margin-top:-1.6rem}}
.stat{{background:#fff;border-radius:14px;padding:1rem;box-shadow:0 8px 22px rgba(16,42,51,.06)}}
.stat b{{display:block;font-family:Fraunces,serif;font-size:1.7rem;color:var(--deep)}}
h2{{font-family:Fraunces,serif}}
.findings{{background:#fff;border-left:4px solid #d4a017;border-radius:12px;padding:1rem 1.2rem}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:1rem}}
.card{{background:#fff;border-radius:14px;padding:.8rem;box-shadow:0 8px 20px rgba(16,42,51,.05)}}
.card img{{width:100%;border-radius:10px;cursor:zoom-in}}
.card h3{{font-size:.92rem;margin:.6rem 0 0}}
.modal{{display:none;position:fixed;inset:0;background:rgba(0,0,0,.88);align-items:center;justify-content:center;z-index:9;padding:1rem;cursor:zoom-out}}
.modal.open{{display:flex}} .modal img{{max-width:95vw;max-height:90vh}}
@media(max-width:700px){{.stats{{grid-template-columns:1fr}}}}
</style></head><body>
<header class="hero"><div>
<p class="brand">Sequence × Embedding</p>
<p>Compare innexin subfamilies with MMseqs2 identity sweeps and protein embeddings (AA composition + 3-mer SVD).</p>
<div class="cta"><a href="#figs">Diagrams</a><a href="../innexin_clade_comparison/index.html">Clade comparison</a><a href="../family_comparison/index.html">vs connexin</a></div>
</div></header>
<main class="wrap">
<div class="stats">
<div class="stat"><b>{n_seq}</b><span>proteins in all-vs-all</span></div>
<div class="stat"><b>{len(THRESHOLDS)}</b><span>MMseqs2 identity thresholds</span></div>
<div class="stat"><b>4</b><span>similarity approaches</span></div>
</div>
<section><h2>Key results</h2><div class="findings"><ol>{finding_li}</ol></div></section>
<section id="figs"><h2>Diagrams</h2><div class="grid">{''.join(cards)}</div></section>
<section><h2>Files</h2>
<p>
<a href="all_innexins.fasta">all_innexins.fasta</a> ·
<a href="sequence_metadata.csv">sequence_metadata.csv</a> ·
<a href="mmseqs_allvsall.m8">mmseqs_allvsall.m8</a> ·
<a href="mmseqs_threshold_sweep_summary.csv">threshold sweep</a> ·
<a href="group_similarity_matrix.csv">group similarity matrix</a> ·
<a href="key_findings.md">key_findings.md</a> ·
<a href="../family_comparison/index.html">innexin vs connexin</a>
</p></section>
</main>
<div id="modal" class="modal" onclick="this.classList.remove('open')"><img id="mimg" alt=""></div>
<script>
function openModal(s){{document.getElementById('mimg').src=s;document.getElementById('modal').classList.add('open')}}
document.addEventListener('keydown',e=>{{if(e.key==='Escape')document.getElementById('modal').classList.remove('open')}})
</script>
</body></html>"""
    (OUT / "index.html").write_text(page, encoding="utf-8")


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    FIGS.mkdir(parents=True, exist_ok=True)
    print("Collecting sequences (refs + discovery + curator translate)…")
    rows = collect_sequences()
    if len(rows) < 5:
        raise SystemExit(f"Too few sequences: {len(rows)}")
    write_fasta_and_meta(rows)
    print(f"Wrote {len(rows)} sequences → {FASTA}")

    print("Running MMseqs2 all-vs-all + cluster sweep…")
    m8, cluster_tsvs = run_mmseqs(rows)
    meta = load_meta_map(rows)
    sweep = summarize_clusters(cluster_tsvs, meta)
    pw = pairwise_from_m8(m8, meta)

    # Embedding approaches
    print("Computing embeddings…")
    aa_mat, ids = aa_composition_matrix(rows)
    aa_sim = cosine_sim(aa_mat)
    km_mat, ids2 = kmer_embedding_matrix(rows, k=3, n_comp=20)
    assert ids == ids2
    km_sim = cosine_sim(km_mat)

    # Group tables
    group_rows: list[dict] = []
    # mmseqs identity group matrix from pairs
    mm_mat_rows = []
    sfs = sorted({r["subfamily"] for r in rows})
    for a in sfs:
        for b in sfs:
            vals = []
            for q, t, pid, sq, st in pw["pairs"]:
                if sq == a and st == b:
                    vals.append(pid)
            if vals:
                mm_mat_rows.append(
                    {
                        "method": "mmseqs_pident",
                        "group_a": a,
                        "group_b": b,
                        "n_pairs": len(vals),
                        "mean_similarity": f"{np.mean(vals):.4f}",
                        "median_similarity": f"{np.median(vals):.4f}",
                        "same_group": str(a == b).lower(),
                    }
                )
    group_rows.extend(mm_mat_rows)
    group_rows.extend(group_similarity_table(aa_sim, ids, meta, "aa_comp_cosine"))
    group_rows.extend(group_similarity_table(km_sim, ids, meta, "kmer3_svd_cosine"))
    write_csv(OUT / "group_similarity_matrix.csv", group_rows)

    # Plots
    plot_threshold_sweep(sweep)
    plot_within_between_identity(pw)
    plot_sf_heatmap(group_rows, "mmseqs_pident", "03_heatmap_mmseqs_identity.png", "Mean MMseqs2 % identity between subfamilies", pct=True)
    plot_sf_heatmap(group_rows, "aa_comp_cosine", "04_heatmap_aa_composition.png", "AA-composition cosine between subfamilies")
    plot_sf_heatmap(group_rows, "kmer3_svd_cosine", "05_heatmap_kmer_embedding.png", "3-mer SVD embedding cosine between subfamilies")
    plot_embedding_pca(km_mat, ids, meta, "3-mer embedding PCA (colored by subfamily)", "06_embedding_pca_kmer.png")
    plot_embedding_pca(aa_mat, ids, meta, "AA-composition PCA (colored by subfamily)", "08_embedding_pca_aa.png")

    emb_within = defaultdict(list)
    id_to_i = {s: i for i, s in enumerate(ids)}
    for sf in sfs:
        members = [r["seq_id"] for r in rows if r["subfamily"] == sf]
        for i, a in enumerate(members):
            for b in members[i + 1 :]:
                emb_within[sf].append(float(km_sim[id_to_i[a], id_to_i[b]]))
    plot_method_agreement(pw["within"], emb_within)

    findings = key_findings(sweep, pw, group_rows)
    (OUT / "key_findings.md").write_text("# Key findings\n\n" + "\n".join(f"- {f}" for f in findings) + "\n", encoding="utf-8")

    fig_plan = [
        ("01_mmseqs_threshold_sweep.png", "MMseqs2 identity-threshold sweep"),
        ("02_within_sf_identity.png", "Within-subfamily sequence identity"),
        ("03_heatmap_mmseqs_identity.png", "Subfamily × subfamily MMseqs2 identity"),
        ("04_heatmap_aa_composition.png", "AA-composition cosine heatmap"),
        ("05_heatmap_kmer_embedding.png", "3-mer embedding cosine heatmap"),
        ("06_embedding_pca_kmer.png", "3-mer embedding PCA"),
        ("07_seq_vs_embedding_cohesion.png", "Sequence vs embedding cohesion"),
        ("08_embedding_pca_aa.png", "AA-composition PCA"),
    ]
    build_html(len(rows), findings, fig_plan)
    print(f"Done → {OUT / 'index.html'}  (n={len(rows)})")


if __name__ == "__main__":
    main()
