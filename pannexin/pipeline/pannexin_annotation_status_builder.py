#!/usr/bin/env python3
"""Classify UniProt/NCBI pannexin annotation status for the starter panel.

Runs tools/classify_family_annotation.py for every organism in
project/metadata/species_config.csv (family=pannexin only), then joins
clade hints and curated-panel copy counts.

Outputs → project/metadata/ + project/results/pannexin_annotation_status/

Does not run SynVoy or genome probes.
"""

from __future__ import annotations

import csv
import html
import re
import subprocess
import sys
from collections import defaultdict
from pathlib import Path

from pipeline.common import METADATA_DIR, PROJECT_ROOT, RESULTS_DIR, read_csv_rows, write_csv

OUT = RESULTS_DIR / "pannexin_annotation_status"
STATUS_CSV = METADATA_DIR / "pannexin_annotation_status.csv"
SUMMARY_CSV = METADATA_DIR / "pannexin_annotation_summary.csv"
CLADE_SUMMARY = OUT / "status_by_clade.csv"
SPECIES_CONFIG = METADATA_DIR / "species_config.csv"
PANEL_META = RESULTS_DIR / "pannexin_panel" / "sequence_metadata.csv"

# Parent thesis tools/ lives one level above the pannexin workspace.
REPO_ROOT = PROJECT_ROOT.parent
CLASSIFIER = REPO_ROOT / "tools" / "classify_family_annotation.py"

CATEGORY_LABEL = {
    "annotated_family_member_found": "annotated",
    "related_or_family_like_entry_found": "related",
    "no_related_entry_found": "none",
}


def _clade_of() -> dict[str, str]:
    out: dict[str, str] = {}
    for row in read_csv_rows(SPECIES_CONFIG):
        org = (row.get("organism") or "").strip()
        if org:
            out[org] = (row.get("clade_hint") or "Other").strip() or "Other"
    return out


def _panel_counts() -> dict[str, dict[str, int]]:
    """organism → {n_proteins, n_panx1, n_panx2, n_panx3, n_other}."""
    counts: dict[str, dict[str, int]] = defaultdict(
        lambda: {"n_proteins": 0, "n_panx1": 0, "n_panx2": 0, "n_panx3": 0, "n_other": 0}
    )
    for row in read_csv_rows(PANEL_META):
        org = (row.get("organism") or "").strip()
        if not org:
            continue
        counts[org]["n_proteins"] += 1
        t = (row.get("panx_type") or "").upper()
        if t == "PANX1":
            counts[org]["n_panx1"] += 1
        elif t == "PANX2":
            counts[org]["n_panx2"] += 1
        elif t == "PANX3":
            counts[org]["n_panx3"] += 1
        else:
            counts[org]["n_other"] += 1
    return dict(counts)


def run_classifier() -> None:
    if not CLASSIFIER.exists():
        raise FileNotFoundError(f"Missing classifier: {CLASSIFIER}")
    METADATA_DIR.mkdir(parents=True, exist_ok=True)
    cmd = [
        sys.executable,
        str(CLASSIFIER),
        "--species-csv",
        str(SPECIES_CONFIG),
        "--families",
        "pannexin",
        "-o",
        str(STATUS_CSV),
        "--summary",
        str(SUMMARY_CSV),
    ]
    print("Running UniProt/NCBI annotation classifier (pannexin panel species)…")
    subprocess.run(cmd, check=True, cwd=str(REPO_ROOT))


FALSE_POSITIVE_HIT = re.compile(
    r"(pknox|pbx|homeobox|knotted|lrrc8|volume[\s_-]?regulated|vrac|gyf protein)",
    re.I,
)


def enrich_and_summarize() -> list[dict[str, str]]:
    clade_of = _clade_of()
    panel = _panel_counts()
    status_rows = read_csv_rows(STATUS_CSV)
    # Index by organism (classifier writes one row per organism×family; we asked pannexin only)
    by_org: dict[str, dict[str, str]] = {}
    for row in status_rows:
        if (row.get("family") or "").strip().lower() != "pannexin":
            continue
        by_org[(row.get("organism") or "").strip()] = row

    enriched: list[dict[str, str]] = []
    for org, clade in sorted(clade_of.items(), key=lambda kv: (kv[1], kv[0])):
        st = by_org.get(org, {})
        pc = panel.get(org, {"n_proteins": 0, "n_panx1": 0, "n_panx2": 0, "n_panx3": 0, "n_other": 0})
        cat = (st.get("final_category") or "").strip()
        term = (st.get("top_matching_term") or "").strip()
        gene = (st.get("top_gene_symbol") or "").strip()
        protein = (st.get("top_protein_name") or "").strip()
        note = ""
        # Classifier sometimes ranks pknox2 / GYF under protein_name:pannexin — filter for display.
        if FALSE_POSITIVE_HIT.search(f"{term} {gene} {protein}"):
            note = f"classifier top hit filtered ({gene or term or protein})"
            if int(pc["n_proteins"]) > 0:
                cat = "annotated_family_member_found"
                term = "panel_curated_pannexin"
            elif cat == "annotated_family_member_found":
                cat = "related_or_family_like_entry_found"
                term = "noisy_hit_filtered"
        elif int(pc["n_proteins"]) > 0 and not cat:
            cat = "annotated_family_member_found"
            term = "panel_curated_pannexin"
        enriched.append(
            {
                "organism": org,
                "clade": clade,
                "final_category": cat or "not_classified",
                "category_short": CATEGORY_LABEL.get(cat, cat or "not_classified"),
                "confidence": (st.get("confidence") or "").strip(),
                "top_matching_term": term,
                "source_db": (st.get("database_source") or st.get("source_db") or "").strip(),
                "n_panel_proteins": str(pc["n_proteins"]),
                "n_panx1": str(pc["n_panx1"]),
                "n_panx2": str(pc["n_panx2"]),
                "n_panx3": str(pc["n_panx3"]),
                "n_other": str(pc["n_other"]),
                "qc_note": note,
                "bottleneck_flag": "yes"
                if clade in ("Jawless vertebrate", "Tunicate / lancelet")
                else "no",
            }
        )

    write_csv(OUT / "species_status.csv", enriched)

    # clade roll-up
    by_clade: dict[str, list[dict[str, str]]] = defaultdict(list)
    for row in enriched:
        by_clade[row["clade"]].append(row)
    clade_rows = []
    for clade, items in sorted(by_clade.items(), key=lambda kv: kv[0]):
        cats = defaultdict(int)
        for r in items:
            cats[r["category_short"]] += 1
        n_prot = sum(int(r["n_panel_proteins"]) for r in items)
        clade_rows.append(
            {
                "clade": clade,
                "n_species_in_config": str(len(items)),
                "n_species_with_panel_protein": str(sum(1 for r in items if int(r["n_panel_proteins"]) > 0)),
                "n_panel_proteins": str(n_prot),
                "median_copies_when_present": _median(
                    [int(r["n_panel_proteins"]) for r in items if int(r["n_panel_proteins"]) > 0]
                ),
                "n_annotated": str(cats.get("annotated", 0)),
                "n_related": str(cats.get("related", 0)),
                "n_none": str(cats.get("none", 0)),
                "species": "; ".join(r["organism"] for r in items),
            }
        )
    write_csv(CLADE_SUMMARY, clade_rows)
    return enriched


def _median(vals: list[int]) -> str:
    if not vals:
        return ""
    s = sorted(vals)
    return str(s[len(s) // 2])


def write_html(rows: list[dict[str, str]]) -> None:
    table = "".join(
        "<tr>"
        f"<td>{html.escape(r['organism'])}</td>"
        f"<td>{html.escape(r['clade'])}</td>"
        f"<td>{html.escape(r['category_short'])}</td>"
        f"<td>{html.escape(r['top_matching_term'])}</td>"
        f"<td>{html.escape(r.get('qc_note') or '')}</td>"
        f"<td>{html.escape(r['n_panel_proteins'])}</td>"
        f"<td>{html.escape(r['n_panx1'])}/{html.escape(r['n_panx2'])}/{html.escape(r['n_panx3'])}</td>"
        f"<td>{html.escape(r['bottleneck_flag'])}</td>"
        "</tr>"
        for r in rows
    )
    clade_rows = read_csv_rows(CLADE_SUMMARY)
    clade_table = "".join(
        "<tr>"
        f"<td>{html.escape(r['clade'])}</td>"
        f"<td>{html.escape(r['n_species_in_config'])}</td>"
        f"<td>{html.escape(r['n_species_with_panel_protein'])}</td>"
        f"<td>{html.escape(r['n_panel_proteins'])}</td>"
        f"<td>{html.escape(r['median_copies_when_present'])}</td>"
        f"<td>{html.escape(r['n_annotated'])}/{html.escape(r['n_related'])}/{html.escape(r['n_none'])}</td>"
        "</tr>"
        for r in clade_rows
    )
    bottlenecks = [r for r in rows if r["bottleneck_flag"] == "yes"]
    bn_html = "".join(
        f"<li><em>{html.escape(r['organism'])}</em> ({html.escape(r['clade'])}): "
        f"{html.escape(r['n_panel_proteins'])} panel protein(s), "
        f"UniProt/NCBI = {html.escape(r['category_short'])}"
        + (f" ({html.escape(r['top_matching_term'])})" if r["top_matching_term"] else "")
        + "</li>"
        for r in bottlenecks
    )
    page = f"""<!DOCTYPE html>
<html lang="en"><head>
<meta charset="utf-8"><meta name="viewport" content="width=device-width, initial-scale=1">
<title>Pannexin annotation status by clade</title>
<link href="https://fonts.googleapis.com/css2?family=Fraunces:opsz,wght@9..144,600&family=Sora:wght@400;600&display=swap" rel="stylesheet">
<style>
:root {{ --ink:#12202a; --muted:#4d6570; --ember:#c2410c; --line:rgba(18,32,42,.12); --bg:#f7f3ef; }}
body {{ margin:0; font-family:Sora,sans-serif; color:var(--ink); background:var(--bg); line-height:1.55 }}
.wrap {{ width:min(1100px, calc(100% - 2rem)); margin:0 auto; padding:2rem 0 4rem }}
h1 {{ font-family:Fraunces,serif; margin:0 0 .4rem; font-size:clamp(1.6rem,3vw,2.2rem) }}
.lede {{ color:var(--muted); max-width:68ch }}
.take {{ border-left:3px solid var(--ember); background:rgba(194,65,12,.07); padding:.75rem .9rem; border-radius:0 12px 12px 0; margin:1rem 0 }}
.table-wrap {{ overflow:auto; background:#fff; border:1px solid var(--line); border-radius:12px; margin:1rem 0 }}
table {{ border-collapse:collapse; width:100%; font-size:.82rem }}
th,td {{ padding:.45rem .55rem; border-bottom:1px solid var(--line); text-align:left }}
th {{ background:rgba(194,65,12,.08) }}
h2 {{ font-family:Fraunces,serif; font-size:1.25rem; margin:1.6rem 0 .5rem }}
.footer a {{ margin-right:.85rem; color:var(--ember) }}
ul {{ color:var(--muted) }}
</style></head><body>
<main class="wrap">
  <h1>Pannexin annotation status</h1>
  <p class="lede">
    Same UniProt/NCBI layered classifier used for innexins and connexins, restricted to the
    pannexin species list. Panel protein counts come from the curated FASTA set.
  </p>
  <p class="take"><strong>Bottleneck check (lancelet / lamprey / tunicate).</strong></p>
  <ul>{bn_html or "<li>No early-chordate species in config.</li>"}</ul>

  <h2>By clade</h2>
  <div class="table-wrap"><table>
    <thead><tr>
      <th>clade</th><th>species</th><th>with panel protein</th><th>proteins</th>
      <th>median copies</th><th>annotated / related / none</th>
    </tr></thead>
    <tbody>{clade_table}</tbody>
  </table></div>

  <h2>Every panel species</h2>
  <div class="table-wrap"><table>
    <thead><tr>
      <th>organism</th><th>clade</th><th>UniProt/NCBI</th><th>term</th><th>QC</th>
      <th>panel n</th><th>1/2/3</th><th>bottleneck?</th>
    </tr></thead>
    <tbody>{table}</tbody>
  </table></div>

  <p class="footer">
    <a href="../pannexin_clade_comparison/index.html">Clade × type</a>
    <a href="../pannexin_curator/index.html">Gaps</a>
    <a href="../pannexin_path/index.html">Path</a>
  </p>
</main></body></html>
"""
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "index.html").write_text(page, encoding="utf-8")


def build(*, skip_classifier: bool = False) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    if skip_classifier and STATUS_CSV.exists():
        print(f"Reusing existing {STATUS_CSV}")
    else:
        run_classifier()
    rows = enrich_and_summarize()
    write_html(rows)
    print(f"Done → {OUT / 'index.html'} ({len(rows)} species)")


if __name__ == "__main__":
    skip = "--reuse" in sys.argv
    build(skip_classifier=skip)
