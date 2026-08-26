#!/usr/bin/env python3
"""Assemble the curated pannexin reference panel from UniProt FASTAs.

Filters false hits (e.g. pknox2), short fragments, and near-duplicate isoforms
(keep longest per organism + gene label). Writes FASTA + metadata tables.

Outputs → project/results/pannexin_panel/
"""

from __future__ import annotations

import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path

from pipeline.common import (
    CLADE_COLORS,
    METADATA_DIR,
    PROJECT_ROOT,
    REFERENCES_DIR,
    RESULTS_DIR,
    TYPE_COLORS,
    classify_pannexin_type,
    read_csv_rows,
    write_csv,
)

OUT = RESULTS_DIR / "pannexin_panel"
FASTA = OUT / "all_pannexins.fasta"
META = OUT / "sequence_metadata.csv"
CLADE_SUMMARY = OUT / "clade_summary.csv"
TYPE_SUMMARY = OUT / "type_summary.csv"

MIN_LEN = 200
FALSE_POSITIVE_GENES = re.compile(
    r"(pknox|pbx|homeobox|knotted|lrrc8|volume[\s_-]?regulated|vrac)",
    re.I,
)


def _clean_aa(seq: str) -> str:
    s = re.sub(r"[^A-Za-z*]", "", seq).upper().split("*")[0]
    return re.sub(r"[^ACDEFGHIKLMNPQRSTVWY]", "X", s)


def _read_one_fasta(path: Path) -> tuple[str, str]:
    hdr, chunks = "", []
    for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
        if line.startswith(">"):
            hdr = line[1:].strip()
        else:
            chunks.append(line.strip())
    return hdr, _clean_aa("".join(chunks))


def _clade_map() -> dict[str, str]:
    out: dict[str, str] = {}
    for row in read_csv_rows(METADATA_DIR / "species_config.csv"):
        org = (row.get("organism") or "").strip()
        if org:
            out[org.replace(" ", "_")] = (row.get("clade_hint") or "Other").strip() or "Other"
            out[org] = out[org.replace(" ", "_")]
    return out


def _gene_label(filename_stem: str, header: str) -> str:
    # filename: Gene__ACCESSION
    gene = filename_stem.split("__")[0] if "__" in filename_stem else filename_stem
    if gene.upper().startswith("LOC") or gene.upper() == "PANX":
        # try header for gene=
        m = re.search(r"GN=(\S+)", header)
        if m:
            gene = m.group(1)
    return gene


def collect_rows() -> list[dict]:
    clade_of = _clade_map()
    rows: list[dict] = []
    if not REFERENCES_DIR.exists():
        return rows

    for path in sorted(REFERENCES_DIR.rglob("*.fasta")):
        organism = path.parent.name
        header, seq = _read_one_fasta(path)
        gene = _gene_label(path.stem, header)
        accession = path.stem.split("__")[-1] if "__" in path.stem else path.stem
        if FALSE_POSITIVE_GENES.search(gene) or FALSE_POSITIVE_GENES.search(header):
            continue
        if len(seq) < MIN_LEN:
            continue
        panx_type = classify_pannexin_type(f"{gene} {header}")
        if panx_type == "other/unknown":
            panx_type = classify_pannexin_type(header)
        clade = clade_of.get(organism, clade_of.get(organism.replace("_", " "), "Other"))
        seq_id = f"panx|{organism}|{gene}|{accession}"
        rows.append(
            {
                "seq_id": seq_id,
                "organism": organism.replace("_", " "),
                "organism_key": organism,
                "gene": gene,
                "accession": accession,
                "panx_type": panx_type,
                "clade": clade,
                "length": str(len(seq)),
                "cys_count": str(seq.count("C")),
                "source_path": str(path.relative_to(PROJECT_ROOT)),
                "seq": seq,
            }
        )
    return rows


def dedupe_longest(rows: list[dict]) -> list[dict]:
    """Keep longest sequence per organism + normalized gene (+ type for generic PANX)."""
    buckets: dict[tuple[str, str], dict] = {}
    for row in rows:
        gene_key = re.sub(r"[^a-z0-9]", "", row["gene"].lower())
        key = (row["organism_key"], gene_key or row["accession"])
        prev = buckets.get(key)
        if prev is None or int(row["length"]) > int(prev["length"]):
            buckets[key] = row
    return sorted(buckets.values(), key=lambda r: (r["clade"], r["organism"], r["panx_type"], r["gene"]))


def write_fasta(rows: list[dict], path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        for row in rows:
            handle.write(f">{row['seq_id']}\n")
            seq = row["seq"]
            for i in range(0, len(seq), 80):
                handle.write(seq[i : i + 80] + "\n")


def summarize(rows: list[dict]) -> None:
    by_clade: dict[str, list[dict]] = defaultdict(list)
    by_type: dict[str, list[dict]] = defaultdict(list)
    for row in rows:
        by_clade[row["clade"]].append(row)
        by_type[row["panx_type"]].append(row)

    clade_rows = []
    for clade, items in sorted(by_clade.items(), key=lambda kv: (-len(kv[1]), kv[0])):
        orgs = sorted({r["organism"] for r in items})
        clade_rows.append(
            {
                "clade": clade,
                "n_proteins": str(len(items)),
                "n_species": str(len(orgs)),
                "species": "; ".join(orgs),
                "color": CLADE_COLORS.get(clade, CLADE_COLORS["Other"]),
            }
        )
    write_csv(CLADE_SUMMARY, clade_rows)

    type_rows = []
    for panx_type, items in sorted(by_type.items()):
        type_rows.append(
            {
                "panx_type": panx_type,
                "n_proteins": str(len(items)),
                "median_length": str(sorted(int(r["length"]) for r in items)[len(items) // 2]),
                "color": TYPE_COLORS.get(panx_type, TYPE_COLORS["other/unknown"]),
            }
        )
    write_csv(TYPE_SUMMARY, type_rows)


def _needs_type_rescue(row: dict) -> bool:
    """Only rescue truly untyped / generic symbols — keep header-typed LOCs."""
    gene = (row.get("gene") or "").upper()
    typ = row.get("panx_type") or "other/unknown"
    if typ == "other/unknown":
        return True
    if gene in {"PANX", "PX"}:
        return True
    return False


def rescue_types_mmseqs(rows: list[dict]) -> list[dict]:
    """Assign PANX1/2/3 to weak/LOC/unknown labels via MMseqs vs typed panel refs."""
    typed = [r for r in rows if r["panx_type"] in {"PANX1", "PANX2", "PANX3"} and not _needs_type_rescue(r)]
    weak = [r for r in rows if _needs_type_rescue(r)]
    for r in rows:
        r.setdefault("panx_type_source", "uniprot_label")
        r.setdefault("type_rescue_hit", "")
        r.setdefault("type_rescue_pident", "")
    if not weak or not typed:
        return rows

    mmseqs = shutil.which("mmseqs") or str(Path.home() / "miniconda3/envs/synvoy_env/bin/mmseqs")
    if not Path(mmseqs).exists():
        print("WARN: mmseqs not found — skip type rescue")
        return rows

    work = OUT / "work_type_rescue"
    if work.exists():
        shutil.rmtree(work)
    work.mkdir(parents=True, exist_ok=True)
    q_fa = work / "query.fasta"
    t_fa = work / "typed.fasta"
    write_fasta(weak, q_fa)
    write_fasta(typed, t_fa)
    m8 = work / "rescue.m8"
    tmp = work / "tmp"
    cmd = [
        mmseqs,
        "easy-search",
        str(q_fa),
        str(t_fa),
        str(m8),
        str(tmp),
        "--threads",
        "1",
        "-e",
        "1e-5",
        "--format-output",
        "query,target,pident,alnlen,evalue",
    ]
    try:
        subprocess.run(cmd, check=True, capture_output=True, text=True)
    except (subprocess.CalledProcessError, FileNotFoundError) as exc:
        print(f"WARN: type rescue MMseqs failed ({exc})")
        return rows

    best: dict[str, tuple[str, float, str]] = {}
    if m8.exists():
        for line in m8.read_text(encoding="utf-8").splitlines():
            parts = line.split("\t")
            if len(parts) < 5:
                continue
            q, t, pident_s, aln_s, _ = parts[:5]
            try:
                pident = float(pident_s)
                aln = float(aln_s)
            except ValueError:
                continue
            if aln < 80:
                continue
            if pident < 40.0:
                continue
            # target type from typed seq_id / meta
            t_type = "other/unknown"
            for tr in typed:
                if tr["seq_id"] == t:
                    t_type = tr["panx_type"]
                    break
            if t_type not in {"PANX1", "PANX2", "PANX3"}:
                continue
            prev = best.get(q)
            if prev is None or pident > prev[1]:
                best[q] = (t_type, pident, t)

    n_rescue = 0
    for r in rows:
        hit = best.get(r["seq_id"])
        if not hit:
            continue
        t_type, pident, t_id = hit
        r["panx_type"] = t_type
        r["panx_type_source"] = "mmseqs_rescue"
        r["type_rescue_hit"] = t_id
        r["type_rescue_pident"] = f"{pident:.1f}"
        n_rescue += 1
    print(f"Type rescue: reassigned {n_rescue}/{len(weak)} weak/LOC/unknown labels via MMseqs")
    return rows


def build() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    raw = collect_rows()
    rows = dedupe_longest(raw)
    rows = rescue_types_mmseqs(rows)
    write_fasta(rows, FASTA)
    meta = [{k: v for k, v in r.items() if k != "seq"} for r in rows]
    write_csv(META, meta)
    summarize(rows)
    print(f"Panel: {len(rows)} proteins from {len({r['organism'] for r in rows})} species → {FASTA}")
    print(f"  (raw kept after length/FP filter: {len(raw)}; after gene dedupe: {len(rows)})")


if __name__ == "__main__":
    build()
