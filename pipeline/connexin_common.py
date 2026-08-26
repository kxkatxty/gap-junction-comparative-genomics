"""Shared connexin typing, clade labels, and sequence I/O."""

from __future__ import annotations

import csv
import re
from collections import Counter
from pathlib import Path

from pipeline.common import METADATA_DIR, PROJECT_ROOT, RESULTS_DIR

REF_ROOT = PROJECT_ROOT / "project/data/references/connexins"
CAND = METADATA_DIR / "gap_junction_candidates_master.csv"
GENE_SUMMARY = RESULTS_DIR / "exon_structures" / "gene_summary.csv"
DISCOVERY_DIRS = (
    RESULTS_DIR / "connexin_evolutionary_discovery",
    RESULTS_DIR / "connexin_evolutionary_batch2_discovery",
    RESULTS_DIR / "connexin_discovery",
)

CLADE_COLORS = {
    "Mammal": "#B45309",
    "Bird": "#0369A1",
    "Reptile": "#15803D",
    "Amphibian": "#0F766E",
    "Ray-finned fish": "#1D4ED8",
    "Cartilaginous fish": "#0E7490",
    "Other vertebrate": "#64748B",
}

TYPE_COLORS = {
    "GJA1": "#C2410C",
    "GJA3": "#EA580C",
    "GJA4": "#D97706",
    "GJA5": "#CA8A04",
    "GJA8": "#A16207",
    "other GJA": "#92400E",
    "GJB1": "#1D4ED8",
    "GJB2": "#2563EB",
    "GJB3": "#0284C7",
    "GJB6": "#0369A1",
    "other GJB": "#075985",
    "GJC": "#15803D",
    "GJD": "#7C3AED",
    "GJE": "#A21CAF",
    "other/unknown": "#94A3B8",
}

SUBFAMILY_COLORS = {
    "alpha_GJA": "#C2410C",
    "beta_GJB": "#1D4ED8",
    "gamma_GJC": "#15803D",
    "delta_GJD": "#7C3AED",
    "epsilon_GJE": "#A21CAF",
    "unassigned": "#94A3B8",
}

SUBFAMILY_FROM_TYPE = {
    "GJA1": "alpha_GJA",
    "GJA3": "alpha_GJA",
    "GJA4": "alpha_GJA",
    "GJA5": "alpha_GJA",
    "GJA8": "alpha_GJA",
    "other GJA": "alpha_GJA",
    "GJB1": "beta_GJB",
    "GJB2": "beta_GJB",
    "GJB3": "beta_GJB",
    "GJB6": "beta_GJB",
    "other GJB": "beta_GJB",
    "GJC": "gamma_GJC",
    "GJD": "delta_GJD",
    "GJE": "epsilon_GJE",
    "other/unknown": "unassigned",
}

BIRD_HINTS = (
    "gallus", "anas", "passer", "luscinia", "cecropis", "butorides", "sturnus",
    "dumetella", "emberiza", "marmaronetta", "aves",
)
MAMMAL_HINTS = (
    "homo", "mus", "rattus", "bos", "sus", "canis", "felis", "equus",
    "oryctolagus", "monodelphis", "ceratotherium", "babyrousa", "antrozous",
    "marmota", "ovis", "pongo", "macaca", "ursus", "hylobates", "mesocricetus",
)
FISH_HINTS = (
    "danio", "oncorhynchus", "gasterosteus", "takifugu", "oreochromis",
    "ictalurus", "hucho", "hemibarbus", "epinephelus", "scleropages",
    "triplophysa", "acanthurus", "actinopteryg",
)
CHONDRICHTHYES_HINTS = (
    "carcharias", "galeorhinus", "echinorhinus", "leucoraja", "malacoraja",
    "okamejei", "potamotrygon", "hydrolagus", "chimaera", "chondrichth",
)
AMPHIBIAN_HINTS = (
    "xenopus", "bombina", "bufotes", "discoglossus", "desmognathus",
    "fejervarya", "gastrophryne", "eleutherodactylus", "epidalea",
    "ichthyophis", "boulenophrys", "amphib",
)
REPTILE_HINTS = (
    "anolis", "python", "chelonia", "correlophus", "gopherus", "mabuya",
    "euprepiophis", "liasis", "caretta", "apalone", "reptil", "squamat",
)


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


def clean_aa(seq: str) -> str:
    s = re.sub(r"[^A-Za-z*]", "", seq).upper().split("*")[0]
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", s)


def read_fasta(path: Path) -> list[tuple[str, str]]:
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


def load_catalog_clades() -> dict[str, str]:
    mapping: dict[str, str] = {}
    coarse = {
        "Mammalia": "Mammal",
        "Mammals": "Mammal",
        "Amphibia": "Amphibian",
        "Amphibians": "Amphibian",
        "Actinopterygii": "Ray-finned fish",
        "Fishes": "Ray-finned fish",
        "Chondrichthyes": "Cartilaginous fish",
        "Birds": "Bird",
        "Reptiles": "Reptile",
        "Sauropsida": "Reptile",
    }
    for path in METADATA_DIR.glob("species_clade_catalog*.csv"):
        with path.open(newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh):
                org = (row.get("organism") or "").strip()
                raw = (row.get("clade") or "").strip()
                notes = (row.get("notes") or "").lower()
                if not org:
                    continue
                clade = coarse.get(raw, "Other vertebrate")
                if raw == "Sauropsida":
                    if "bird" in notes or any(h in org.lower() for h in BIRD_HINTS):
                        clade = "Bird"
                    else:
                        clade = "Reptile"
                mapping[org] = clade
                mapping[org.replace(" ", "_")] = clade
    return mapping


_CATALOG = None


def catalog() -> dict[str, str]:
    global _CATALOG
    if _CATALOG is None:
        _CATALOG = load_catalog_clades()
    return _CATALOG


def infer_clade(species_dir: str, organism: str = "") -> str:
    cat = catalog()
    if organism in cat:
        return cat[organism]
    if species_dir in cat:
        return cat[species_dir]
    name = (organism or species_dir or "").lower().replace("_", " ")
    if any(h in name for h in MAMMAL_HINTS):
        return "Mammal"
    if any(h in name for h in BIRD_HINTS):
        return "Bird"
    if any(h in name for h in AMPHIBIAN_HINTS):
        return "Amphibian"
    if any(h in name for h in REPTILE_HINTS):
        return "Reptile"
    if any(h in name for h in CHONDRICHTHYES_HINTS):
        return "Cartilaginous fish"
    if any(h in name for h in FISH_HINTS):
        return "Ray-finned fish"
    return "Other vertebrate"


def classify_connexin_type(text: str) -> str:
    t = (text or "").upper().replace("-", "").replace("_", " ")
    rules = [
        (r"GJA10|CX62", "other GJA"),
        (r"GJA9|CX59", "other GJA"),
        (r"GJA8|CX50", "GJA8"),
        (r"GJA6", "other GJA"),
        (r"GJA5|CX40", "GJA5"),
        (r"GJA4|CX37", "GJA4"),
        (r"GJA3|CX46", "GJA3"),
        (r"GJA1|CX43", "GJA1"),
        (r"GJA", "other GJA"),
        (r"GJB7", "other GJB"),
        (r"GJB6|CX30", "GJB6"),
        (r"GJB5", "other GJB"),
        (r"GJB4", "other GJB"),
        (r"GJB3|CX31", "GJB3"),
        (r"GJB2|CX26", "GJB2"),
        (r"GJB1|CX32", "GJB1"),
        (r"GJB", "other GJB"),
        (r"GJC", "GJC"),
        (r"GJD", "GJD"),
        (r"GJE", "GJE"),
    ]
    for pat, label in rules:
        if re.search(pat, t):
            return label
    return "other/unknown"


def is_kept_rank(rank: str) -> bool:
    r = (rank or "").lower()
    if "reject" in r:
        return False
    return any(
        k in r
        for k in (
            "high_confidence",
            "accepted",
            "weak_manual_review",
            "probable",
            "present",
        )
    )


def is_high_conf(rank: str) -> bool:
    r = (rank or "").lower()
    return "high_confidence" in r or r in {"accepted", "reference"}


def load_reference_rows() -> list[dict]:
    rows: list[dict] = []
    if not REF_ROOT.exists():
        return rows
    for fasta in sorted(REF_ROOT.rglob("*.fasta")):
        organism = fasta.parent.name.replace("_", " ")
        species_dir = fasta.parent.name
        gene = fasta.stem.split("__")[0]
        acc = fasta.stem.split("__")[1] if "__" in fasta.stem else ""
        seq = "".join(s for _, s in read_fasta(fasta))
        aa = clean_aa(seq)
        rtype = classify_connexin_type(gene)
        rows.append(
            {
                "source": "reference_db",
                "species_dir": species_dir,
                "organism": organism,
                "clade": infer_clade(species_dir, organism),
                "candidate_id": fasta.stem,
                "gene_label": gene,
                "reference_type": rtype,
                "subfamily": SUBFAMILY_FROM_TYPE.get(rtype, "unassigned"),
                "protein_length": len(aa),
                "exon_count": "",
                "tm_helix_count": "",
                "reference_identity": 1.0,
                "best_reference_hit": f"{gene}|{acc}" if acc else gene,
                "rank_category": "reference",
                "discovery_batch": "reference_db",
                "seq": aa,
                "seq_id": f"ref__{species_dir}__{gene}__{acc or 'na'}",
            }
        )
    return rows


def load_candidate_table_rows() -> list[dict]:
    rows: list[dict] = []
    if not CAND.exists():
        return rows
    with CAND.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("family") != "connexin":
                continue
            rank = r.get("rank_category") or ""
            if not is_kept_rank(rank):
                continue
            species_dir = r.get("species_dir") or ""
            organism = r.get("organism") or species_dir.replace("_", " ")
            hit = r.get("best_reference_hit") or ""
            rtype = classify_connexin_type(hit or r.get("candidate_id") or "")
            batch = r.get("discovery_batch") or ""
            source = "discovery" if "discover" in batch else "discovery"
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
                    "seq": "",
                    "seq_id": f"disc__{species_dir}__{r.get('candidate_id') or 'na'}",
                }
            )
    return rows


def load_discovery_fastas(high_conf_only: bool = True) -> list[dict]:
    """Parse candidate_connexins.fasta files and keep plausible products."""
    rows: list[dict] = []
    seen: set[str] = set()
    for root in DISCOVERY_DIRS:
        if not root.exists():
            continue
        for fasta in sorted(root.glob("*/candidate_connexins.fasta")):
            species_dir = fasta.parent.name
            organism = species_dir.replace("_", " ")
            for hdr, seq in read_fasta(fasta):
                aa = clean_aa(seq)
                rank = ""
                m = re.search(r"rank=(\S+)", hdr)
                if m:
                    rank = m.group(1)
                rlow = rank.lower()
                if high_conf_only:
                    keep = is_high_conf(rank) or "possible_connexin" in rlow
                    if not keep:
                        continue
                    if is_high_conf(rank) and not (150 <= len(aa) <= 520):
                        continue
                    if "possible_connexin" in rlow and not (200 <= len(aa) <= 450):
                        continue
                else:
                    if not is_kept_rank(rank):
                        continue
                    if not (150 <= len(aa) <= 520):
                        continue
                cid = hdr.split()[0]
                key = f"{species_dir}|{cid}|{len(aa)}"
                if key in seen:
                    continue
                seen.add(key)
                rtype = classify_connexin_type(hdr)
                ex = re.search(r"exons=(\d+)", hdr)
                tm = re.search(r"tm=(\d+)", hdr)
                rows.append(
                    {
                        "source": "discovery",
                        "species_dir": species_dir,
                        "organism": organism,
                        "clade": infer_clade(species_dir, organism),
                        "candidate_id": cid,
                        "gene_label": cid,
                        "reference_type": rtype,
                        "subfamily": SUBFAMILY_FROM_TYPE.get(rtype, "unassigned"),
                        "protein_length": len(aa),
                        "exon_count": _i(ex.group(1) if ex else 0),
                        "tm_helix_count": _i(tm.group(1) if tm else 0),
                        "reference_identity": 0.0,
                        "best_reference_hit": "",
                        "rank_category": rank or "discovery",
                        "discovery_batch": fasta.parent.parent.name,
                        "seq": aa,
                        "seq_id": f"disc__{species_dir}__{cid}",
                    }
                )
    return rows


def load_exon_rows() -> list[dict]:
    rows: list[dict] = []
    if not GENE_SUMMARY.exists():
        return rows
    with GENE_SUMMARY.open(newline="", encoding="utf-8") as fh:
        for r in csv.DictReader(fh):
            if r.get("family") != "connexin":
                continue
            gene = r.get("reference_gene_symbol") or r.get("gene_symbol") or ""
            organism = r.get("organism") or ""
            species_dir = organism.replace(" ", "_")
            rtype = classify_connexin_type(gene)
            rows.append(
                {
                    "organism": organism,
                    "species_dir": species_dir,
                    "clade": infer_clade(species_dir, organism),
                    "gene_label": gene,
                    "reference_type": rtype,
                    "subfamily": SUBFAMILY_FROM_TYPE.get(rtype, "unassigned"),
                    "mean_exon_count": _f(r.get("mean_exon_count")),
                    "gene_span": _f(r.get("gene_span")),
                }
            )
    return rows


def write_fasta(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as fh:
        for r in rows:
            seq = clean_aa(r.get("seq") or "")
            if len(seq) < 80:
                continue
            fh.write(f">{r['seq_id']}\n")
            for i in range(0, len(seq), 80):
                fh.write(seq[i : i + 80] + "\n")


def apply_best_hits(rows: list[dict], m8_path: Path) -> None:
    """Overwrite discovery types from MMseqs best hit vs reference connexins."""
    if not m8_path.exists():
        return
    best: dict[str, tuple[float, str, float]] = {}
    for line in m8_path.read_text(encoding="utf-8").splitlines():
        parts = line.split("\t")
        if len(parts) < 3:
            continue
        q, t, pident = parts[0], parts[1], float(parts[2])
        prev = best.get(q)
        bits = float(parts[11]) if len(parts) > 11 else pident
        if prev is None or bits > prev[0]:
            best[q] = (bits, t, pident)
    for r in rows:
        hit = best.get(r["seq_id"])
        if not hit:
            continue
        _, target, pident = hit
        r["best_reference_hit"] = target
        r["reference_identity"] = pident / 100.0 if pident > 1.5 else pident
        rtype = classify_connexin_type(target)
        r["reference_type"] = rtype
        r["subfamily"] = SUBFAMILY_FROM_TYPE.get(rtype, "unassigned")
        # gene label from target id like ref__Homo_sapiens__GJA1__P17302
        bits = target.split("__")
        if len(bits) >= 3:
            r["gene_label"] = bits[2]


def summarize(rows: list[dict]) -> dict:
    return {
        "n_loci": len(rows),
        "n_species": len({r["species_dir"] for r in rows}),
        "by_clade": Counter(r["clade"] for r in rows),
        "by_type": Counter(r["reference_type"] for r in rows),
        "by_source": Counter(r["source"] for r in rows),
        "by_sf": Counter(r["subfamily"] for r in rows),
    }
