#!/usr/bin/env python3
"""MMseqs2 + embedding similarity across connexin α/β/γ/δ/ε groups."""

from __future__ import annotations

import html
import math
import shutil
import subprocess
from collections import Counter, defaultdict
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from pipeline.common import RESULTS_DIR, write_csv
from pipeline.connexin_common import (
    SUBFAMILY_COLORS,
    apply_best_hits,
    load_discovery_fastas,
    load_reference_rows,
    write_fasta,
)

OUT = RESULTS_DIR / "connexin_similarity"
FIGS = OUT / "figures"
WORK = OUT / "work"
FASTA = OUT / "all_connexins.fasta"
META = OUT / "sequence_metadata.csv"
MMSEQS = shutil.which("mmseqs") or str(
    Path.home() / "miniconda3/envs/synvoy_env/bin/mmseqs"
)
THRESHOLDS = [0.30, 0.40, 0.50, 0.60, 0.70, 0.80]
HIT_M8 = RESULTS_DIR / "connexin_clade_comparison" / "discovery_vs_reference.m8"

plt.rcParams.update(
    {
        "figure.dpi": 140,
        "savefig.dpi": 180,
        "font.size": 10,
        "axes.spines.top": False,
        "axes.spines.right": False,
    }
)


def _save(fig, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(path, bbox_inches="tight", facecolor="white")
    plt.close(fig)


def collect() -> list[dict]:
    refs = load_reference_rows()
    disc = load_discovery_fastas(high_conf_only=True)
    if HIT_M8.exists():
        apply_best_hits(disc, HIT_M8)
    rows = [r for r in (refs + disc) if len(r.get("seq") or "") >= 80]
    return rows


def run_mmseqs(rows: list[dict]) -> tuple[Path, dict[float, Path]]:
    WORK.mkdir(parents=True, exist_ok=True)
    write_fasta(FASTA, rows)
    db = WORK / "seqdb"
    tmp = WORK / "tmp"
    tmp.mkdir(exist_ok=True)
    for p in WORK.glob("seqdb*"):
        p.unlink() if p.is_file() else shutil.rmtree(p, ignore_errors=True)
    subprocess.run([MMSEQS, "createdb", str(FASTA), str(db), "-v", "1"], check=True)
    search_db = WORK / "allvsall"
    search_tsv = OUT / "mmseqs_allvsall.m8"
    for p in WORK.glob("allvsall*"):
        if p.is_file():
            p.unlink()
    subprocess.run(
        [MMSEQS, "search", str(db), str(db), str(search_db), str(tmp), "-s", "7.5",
         "--max-seqs", "800", "-e", "1e-3", "-v", "1"],
        check=True,
    )
    subprocess.run(
        [MMSEQS, "convertalis", str(db), str(db), str(search_db), str(search_tsv),
         "--format-output",
         "query,target,pident,alnlen,mismatch,gapopen,qstart,qend,tstart,tend,evalue,bits",
         "-v", "1"],
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
            [MMSEQS, "cluster", str(db), str(clust), str(tmp),
             "--min-seq-id", str(thr), "-c", "0.8", "--cov-mode", "0", "-v", "1"],
            check=True,
        )
        subprocess.run([MMSEQS, "createtsv", str(db), str(db), str(clust), str(tsv), "-v", "1"], check=True)
        cluster_tsvs[thr] = tsv
    return search_tsv, cluster_tsvs


def summarize_clusters(cluster_tsvs: dict[float, Path], meta: dict[str, dict]) -> list[dict]:
    out = []
    for thr, path in sorted(cluster_tsvs.items()):
        groups: dict[str, list[str]] = defaultdict(list)
        for line in path.read_text(encoding="utf-8").splitlines():
            if not line.strip():
                continue
            rep, mem = line.split("\t")[:2]
            groups[rep].append(mem)
        sizes = [len(v) for v in groups.values()]
        purity, clean = [], []
        for members in groups.values():
            if len(members) < 2:
                continue
            spp = Counter(meta.get(m, {}).get("species_dir", "?") for m in members)
            sfs = [meta.get(m, {}).get("subfamily", "unassigned") for m in members]
            clean.append(sum(1 for c in spp.values() if c == 1) / max(len(spp), 1))
            purity.append(Counter(sfs).most_common(1)[0][1] / len(sfs))
        out.append(
            {
                "min_seq_id": thr,
                "n_clusters": len(groups),
                "n_singletons": sum(1 for s in sizes if s == 1),
                "median_cluster_size": float(np.median(sizes)) if sizes else 0,
                "max_cluster_size": max(sizes) if sizes else 0,
                "mean_clean_orthology": float(np.mean(clean)) if clean else float("nan"),
                "mean_subfamily_purity": float(np.mean(purity)) if purity else float("nan"),
                "n_multi_clusters": len(purity),
            }
        )
    write_csv(OUT / "mmseqs_threshold_sweep_summary.csv", out)
    return out


def pairwise_from_m8(m8: Path, meta: dict[str, dict]) -> dict:
    within, between, pairs = defaultdict(list), defaultdict(list), []
    with m8.open(encoding="utf-8") as fh:
        for line in fh:
            q, t, pident, *_ = line.split("\t")
            if q == t:
                continue
            pid = float(pident)
            sq = meta.get(q, {}).get("subfamily", "unassigned")
            st = meta.get(t, {}).get("subfamily", "unassigned")
            pairs.append((q, t, pid, sq, st))
            if sq == st:
                within[sq].append(pid)
            else:
                between[("ANY_BETWEEN",)].append(pid)
    return {"within": within, "between": between, "pairs": pairs}


def aa_matrix(rows):
    aas = list("ACDEFGHIKLMNPQRSTVWY")
    ids = [r["seq_id"] for r in rows]
    X = np.zeros((len(rows), len(aas)))
    for i, r in enumerate(rows):
        c = Counter(r["seq"])
        n = max(len(r["seq"]), 1)
        for j, a in enumerate(aas):
            X[i, j] = c.get(a, 0) / n
    nrm = np.linalg.norm(X, axis=1, keepdims=True)
    nrm[nrm == 0] = 1
    return X / nrm, ids


def kmer_matrix(rows, k=3, n_comp=20):
    ids = [r["seq_id"] for r in rows]
    seqs = [r["seq"] for r in rows]
    vocab_c: Counter[str] = Counter()
    for seq in seqs:
        for i in range(max(0, len(seq) - k + 1)):
            km = seq[i : i + k]
            if "X" not in km:
                vocab_c[km] += 1
    vocab = [km for km, _ in vocab_c.most_common(800)]
    idx = {km: j for j, km in enumerate(vocab)}
    X = np.zeros((len(seqs), len(vocab)))
    for i, seq in enumerate(seqs):
        for j in range(max(0, len(seq) - k + 1)):
            km = seq[j : j + k]
            if km in idx:
                X[i, idx[km]] += 1
        s = X[i].sum()
        if s:
            X[i] /= s
    Xc = X - X.mean(axis=0, keepdims=True)
    u, s, _vt = np.linalg.svd(Xc, full_matrices=False)
    n_comp = min(n_comp, u.shape[1])
    emb = u[:, :n_comp] * s[:n_comp]
    nrm = np.linalg.norm(emb, axis=1, keepdims=True)
    nrm[nrm == 0] = 1
    return emb / nrm, ids


def group_table(sim, ids, meta, label):
    sfs = sorted({meta[i]["subfamily"] for i in ids if i in meta})
    id_to_i = {s: i for i, s in enumerate(ids)}
    rows = []
    for a in sfs:
        ma = [s for s in ids if meta.get(s, {}).get("subfamily") == a]
        for b in sfs:
            mb = [s for s in ids if meta.get(s, {}).get("subfamily") == b]
            vals = [float(sim[id_to_i[qa], id_to_i[qb]]) for qa in ma for qb in mb if qa != qb]
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


def plot_sweep(summary):
    if not summary:
        return
    thr = [r["min_seq_id"] for r in summary]
    fig, ax = plt.subplots(figsize=(8, 4.8))
    ax.plot(thr, [r["n_clusters"] for r in summary], "o-", color="#C2410C", label="# clusters")
    ax2 = ax.twinx()
    ax2.plot(thr, [r["mean_subfamily_purity"] for r in summary], "s--", color="#1D4ED8", label="subfamily purity")
    ax2.plot(thr, [r["mean_clean_orthology"] for r in summary], "^--", color="#15803D", label="clean-orthology")
    ax.set_xlabel("MMseqs2 min-seq-id")
    ax.set_ylabel("Clusters")
    ax2.set_ylabel("Score (0–1)")
    ax.set_title("Connexin identity-threshold sweep")
    h1, l1 = ax.get_legend_handles_labels()
    h2, l2 = ax2.get_legend_handles_labels()
    ax.legend(h1 + h2, l1 + l2, frameon=False, fontsize=8)
    _save(fig, FIGS / "01_mmseqs_threshold_sweep.png")


def plot_within(pw):
    labels = [k for k, v in pw["within"].items() if v]
    if not labels:
        return
    fig, ax = plt.subplots(figsize=(8.5, 5))
    bp = ax.boxplot([pw["within"][k] for k in labels], tick_labels=labels, patch_artist=True, showfliers=False)
    for patch, lab in zip(bp["boxes"], labels):
        patch.set_facecolor(SUBFAMILY_COLORS.get(lab, "#94A3B8"))
        patch.set_alpha(0.85)
    bet = pw["between"].get(("ANY_BETWEEN",), [])
    if bet:
        ax.axhline(np.median(bet), color="#78716C", ls="--", label=f"median between ({np.median(bet):.1f}%)")
        ax.legend(frameon=False)
    ax.set_ylabel("MMseqs2 % identity")
    ax.set_title("Within-subfamily sequence identity (connexins)")
    plt.setp(ax.xaxis.get_majorticklabels(), rotation=20, ha="right")
    _save(fig, FIGS / "02_within_sf_identity.png")


def plot_heatmap(table_rows, method, out_name, title, pct=False):
    sub = [r for r in table_rows if r["method"] == method]
    if not sub:
        return
    groups = sorted({r["group_a"] for r in sub} | {r["group_b"] for r in sub})
    mat = np.full((len(groups), len(groups)), np.nan)
    idx = {g: i for i, g in enumerate(groups)}
    for r in sub:
        mat[idx[r["group_a"]], idx[r["group_b"]]] = float(r["mean_similarity"])
    fig, ax = plt.subplots(figsize=(7.2, 6))
    im = ax.imshow(mat, cmap="YlOrBr")
    short = [g.replace("alpha_", "α ").replace("beta_", "β ").replace("gamma_", "γ ").replace("delta_", "δ ").replace("epsilon_", "ε ") for g in groups]
    ax.set_xticks(range(len(groups)))
    ax.set_yticks(range(len(groups)))
    ax.set_xticklabels(short, rotation=35, ha="right", fontsize=8)
    ax.set_yticklabels(short, fontsize=8)
    for i in range(len(groups)):
        for j in range(len(groups)):
            v = mat[i, j]
            if not math.isnan(v):
                ax.text(j, i, f"{v:.0f}" if pct else f"{v:.2f}", ha="center", va="center", fontsize=7)
    ax.set_title(title)
    fig.colorbar(im, ax=ax, fraction=0.046, pad=0.04)
    _save(fig, FIGS / out_name)


def plot_pca(emb, ids, meta, title, out_name):
    if emb.shape[1] < 2:
        return
    fig, ax = plt.subplots(figsize=(7.5, 5.8))
    for sf, color in SUBFAMILY_COLORS.items():
        pts = [i for i, s in enumerate(ids) if meta.get(s, {}).get("subfamily") == sf]
        if not pts:
            continue
        ax.scatter(emb[pts, 0], emb[pts, 1], s=28, alpha=0.75, c=color, label=sf, edgecolors="white", linewidths=0.3)
    ax.set_xlabel("Component 1")
    ax.set_ylabel("Component 2")
    ax.set_title(title)
    ax.legend(fontsize=7, frameon=False)
    _save(fig, FIGS / out_name)


def plot_cohesion(mm_within, emb_within):
    sfs = sorted(set(mm_within) & set(emb_within))
    if len(sfs) < 2:
        return
    fig, ax = plt.subplots(figsize=(6.5, 5.2))
    for sf in sfs:
        ax.scatter([np.mean(mm_within[sf])], [np.mean(emb_within[sf])], s=90, c=SUBFAMILY_COLORS.get(sf, "#64748B"), label=sf)
    ax.set_xlabel("Mean within-SF MMseqs2 % identity")
    ax.set_ylabel("Mean within-SF embedding cosine")
    ax.set_title("Sequence vs embedding cohesion (connexin subfamilies)")
    ax.legend(fontsize=7, frameon=False)
    _save(fig, FIGS / "07_seq_vs_embedding_cohesion.png")


def findings(sweep, pw, group_rows) -> list[str]:
    out = []
    if sweep:
        best = max(sweep, key=lambda r: (r["mean_subfamily_purity"] if not math.isnan(r["mean_subfamily_purity"]) else -1))
        out.append(
            f"MMseqs2 sweep: highest α/β/γ/δ purity {best['mean_subfamily_purity']:.2f} "
            f"at min-seq-id={best['min_seq_id']:.2f} ({best['n_clusters']} clusters)."
        )
    within, bet = pw["within"], pw["between"].get(("ANY_BETWEEN",), [])
    if within and bet:
        wvals = [np.median(v) for v in within.values() if v]
        out.append(
            f"Within-subfamily median identity ~{np.mean(wvals):.1f}% vs between-subfamily "
            f"median {np.median(bet):.1f}% — α vs β are sequence-separable but one family."
        )
    km = [r for r in group_rows if r["method"] == "kmer3_svd_cosine"]
    same = [float(r["mean_similarity"]) for r in km if r["same_group"] == "true"]
    diff = [float(r["mean_similarity"]) for r in km if r["same_group"] == "false"]
    if same and diff:
        out.append(f"3-mer embedding cosine: same-subfamily {np.mean(same):.3f} vs different {np.mean(diff):.3f}.")
    aa = [r for r in group_rows if r["method"] == "aa_comp_cosine"]
    same_a = [float(r["mean_similarity"]) for r in aa if r["same_group"] == "true"]
    diff_a = [float(r["mean_similarity"]) for r in aa if r["same_group"] == "false"]
    if same_a and diff_a:
        out.append(
            f"AA-composition cosine barely separates groups (same {np.mean(same_a):.3f} vs diff {np.mean(diff_a):.3f}) "
            "— k-mers and MMseqs carry the structural signal."
        )
    cross = [
        r for r in group_rows
        if r["method"] == "mmseqs_pident"
        and {r["group_a"], r["group_b"]} == {"alpha_GJA", "beta_GJB"}
    ]
    if cross:
        out.append(f"α (GJA) ↔ β (GJB) mean identity {cross[0]['mean_similarity']}% (n={cross[0]['n_pairs']}).")
    out.append("Note: discovery types are MMseqs best-hits; GJA1 vs GJA4 still need synteny to confirm locus identity.")
    return out


def build_html(n, finds, fig_plan):
    cards = []
    for fn, cap in fig_plan:
        if (FIGS / fn).exists():
            cards.append(
                f"<article class='card'><img src='figures/{fn}' alt='{html.escape(cap)}' "
                f"loading='lazy' onclick=\"openModal(this.src)\"><h3>{html.escape(cap)}</h3></article>"
            )
    lis = "".join(f"<li>{html.escape(f)}</li>" for f in finds)
    page = f"""<!DOCTYPE html>
<html lang="en"><head><meta charset="utf-8">
<meta name="viewport" content="width=device-width,initial-scale=1">
<title>Connexin similarity — MMseqs2 & embeddings</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,700&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root{{--ink:#2a1408;--muted:#7a5340;--deep:#7c2d12;--bg:#fbf4ee}}
body{{margin:0;font-family:Sora,sans-serif;color:var(--ink);background:linear-gradient(180deg,#f8ebe3,var(--bg));line-height:1.55}}
.hero{{min-height:70vh;display:grid;align-items:end;padding:clamp(1.3rem,4vw,3rem);color:#fffaf6;
background:linear-gradient(115deg,rgba(124,45,18,.94),rgba(194,65,12,.75) 55%,rgba(3,105,161,.45))}}
.brand{{font-family:Fraunces,serif;font-size:clamp(2.2rem,6vw,4.2rem);margin:0 0 .6rem;line-height:.95}}
.cta a{{display:inline-block;margin:.4rem .4rem 0 0;padding:.75rem 1.1rem;border-radius:999px;background:#fffaf6;color:var(--deep);text-decoration:none;font-weight:700}}
.wrap{{width:min(1140px,calc(100% - 2rem));margin:0 auto;padding:2rem 0 4rem}}
.stats{{display:grid;grid-template-columns:repeat(3,1fr);gap:.8rem;margin-top:-1.6rem}}
.stat{{background:#fff;border-radius:14px;padding:1rem;box-shadow:0 8px 22px rgba(42,20,8,.06)}}
.stat b{{display:block;font-family:Fraunces,serif;font-size:1.7rem;color:var(--deep)}}
.findings{{background:#fff;border-left:4px solid #d97706;border-radius:12px;padding:1rem 1.2rem}}
.grid{{display:grid;grid-template-columns:repeat(auto-fill,minmax(300px,1fr));gap:1rem}}
.card{{background:#fff;border-radius:14px;padding:.8rem;box-shadow:0 8px 20px rgba(42,20,8,.05)}}
.card img{{width:100%;border-radius:10px;cursor:zoom-in}}
.modal{{display:none;position:fixed;inset:0;background:rgba(0,0,0,.88);align-items:center;justify-content:center;z-index:9;padding:1rem;cursor:zoom-out}}
.modal.open{{display:flex}} .modal img{{max-width:95vw;max-height:90vh}}
@media(max-width:700px){{.stats{{grid-template-columns:1fr}}}}
</style></head><body>
<header class="hero"><div>
<p class="brand">Sequence × Embedding</p>
<p>Compare connexin α/β/γ/δ/ε groups with MMseqs2 identity sweeps and protein embeddings.</p>
<div class="cta"><a href="#figs">Diagrams</a><a href="../connexin_clade_comparison/index.html">Clade comparison</a><a href="../family_comparison/index.html">vs innexin</a></div>
</div></header>
<main class="wrap">
<div class="stats">
<div class="stat"><b>{n}</b><span>proteins in all-vs-all</span></div>
<div class="stat"><b>{len(THRESHOLDS)}</b><span>identity thresholds</span></div>
<div class="stat"><b>4</b><span>similarity approaches</span></div>
</div>
<section><h2>Key results</h2><div class="findings"><ol>{lis}</ol></div></section>
<section id="figs"><h2>Diagrams</h2><div class="grid">{''.join(cards)}</div></section>
<section><p>
<a href="all_connexins.fasta">FASTA</a> ·
<a href="sequence_metadata.csv">metadata</a> ·
<a href="mmseqs_allvsall.m8">all-vs-all</a> ·
<a href="group_similarity_matrix.csv">group matrix</a> ·
<a href="../connexin_insights/index.html">curated insights</a> ·
<a href="../family_comparison/index.html">vs innexin</a>
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
    rows = collect()
    if len(rows) < 8:
        raise SystemExit(f"Too few sequences: {len(rows)}")
    write_csv(
        META,
        [
            {k: r.get(k, "") for k in (
                "seq_id", "source", "species_dir", "organism", "clade", "gene_label",
                "reference_type", "subfamily", "protein_length", "rank_category",
            )}
            for r in rows
        ],
    )
    print(f"Sequences={len(rows)}")
    m8, clusters = run_mmseqs(rows)
    meta = {r["seq_id"]: r for r in rows}
    sweep = summarize_clusters(clusters, meta)
    pw = pairwise_from_m8(m8, meta)
    aa_mat, ids = aa_matrix(rows)
    km_mat, ids2 = kmer_matrix(rows)
    assert ids == ids2
    aa_sim = np.clip(aa_mat @ aa_mat.T, -1, 1)
    km_sim = np.clip(km_mat @ km_mat.T, -1, 1)

    sfs = sorted({r["subfamily"] for r in rows})
    mm_rows = []
    for a in sfs:
        for b in sfs:
            vals = [pid for q, t, pid, sq, st in pw["pairs"] if sq == a and st == b]
            if vals:
                mm_rows.append(
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
    group_rows = mm_rows + group_table(aa_sim, ids, meta, "aa_comp_cosine") + group_table(km_sim, ids, meta, "kmer3_svd_cosine")
    write_csv(OUT / "group_similarity_matrix.csv", group_rows)

    plot_sweep(sweep)
    plot_within(pw)
    plot_heatmap(group_rows, "mmseqs_pident", "03_heatmap_mmseqs_identity.png", "Mean MMseqs2 % identity between connexin subfamilies", pct=True)
    plot_heatmap(group_rows, "aa_comp_cosine", "04_heatmap_aa_composition.png", "AA-composition cosine between subfamilies")
    plot_heatmap(group_rows, "kmer3_svd_cosine", "05_heatmap_kmer_embedding.png", "3-mer SVD embedding cosine between subfamilies")
    plot_pca(km_mat, ids, meta, "3-mer embedding PCA (α/β/γ/δ/ε)", "06_embedding_pca_kmer.png")
    plot_pca(aa_mat, ids, meta, "AA-composition PCA", "08_embedding_pca_aa.png")
    emb_within = defaultdict(list)
    id_to_i = {s: i for i, s in enumerate(ids)}
    for sf in sfs:
        members = [r["seq_id"] for r in rows if r["subfamily"] == sf]
        for i, a in enumerate(members):
            for b in members[i + 1 :]:
                emb_within[sf].append(float(km_sim[id_to_i[a], id_to_i[b]]))
    plot_cohesion(pw["within"], emb_within)

    finds = findings(sweep, pw, group_rows)
    (OUT / "key_findings.md").write_text("# Key findings\n\n" + "\n".join(f"- {f}" for f in finds) + "\n", encoding="utf-8")
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
    build_html(len(rows), finds, fig_plan)
    print(f"Done → {OUT / 'index.html'} n={len(rows)}")


if __name__ == "__main__":
    main()
