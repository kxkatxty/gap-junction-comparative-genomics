from __future__ import annotations

import argparse
import csv
import html
import json
import re
import sys
import time
from collections import Counter
from pathlib import Path
from typing import Optional
from urllib.parse import quote
from urllib.request import Request, urlopen

USER_AGENT = "gap-junction-quality-table/1.0"

PROJECT_ROOT = Path("project")
REFERENCES_DIR = PROJECT_ROOT / "data" / "references"
OUTPUT_DIR = PROJECT_ROOT / "metadata"
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)


def http_get_json(url: str) -> dict:
    req = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(req, timeout=60) as resp:
        return json.loads(resp.read().decode("utf-8"))


def find_fasta_files() -> list[tuple[str, Path]]:
    files: list[tuple[str, Path]] = []
    for family in ["innexins", "connexins", "pannexins"]:
        folder = REFERENCES_DIR / family
        if folder.exists():
            for path in sorted(folder.rglob("*.fasta")):
                files.append((family[:-1], path))  # innexins -> innexin
    return files


def read_first_fasta_header(path: Path) -> str:
    with path.open("r", encoding="utf-8") as f:
        for line in f:
            if line.startswith(">"):
                return line.strip()
    return ""


def extract_uniprot_accession(path: Path, header: str) -> Optional[str]:
    match = re.search(r"__([A-Z0-9]{6,10})\.fasta$", path.name)
    if match:
        return match.group(1)

    match = re.search(r"^>\w+\|([A-Z0-9]{6,10})\|", header)
    if match:
        return match.group(1)

    match = re.search(r"\b([A-Z0-9]{6,10})\b", header)
    if match:
        return match.group(1)

    return None


def get_uniprot_entry(accession: str) -> dict:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.json"
    return http_get_json(url)


def parse_uniprot(entry: dict) -> dict:
    accession = entry.get("primaryAccession", "")
    organism = entry.get("organism", {}).get("scientificName", "")
    entry_type = entry.get("entryType", "")
    reviewed = "reviewed" if "reviewed" in entry_type.lower() else "unreviewed"

    protein_desc = entry.get("proteinDescription", {})
    protein_name = ""
    rec = protein_desc.get("recommendedName")
    if rec:
        protein_name = rec.get("fullName", {}).get("value", "")
    if not protein_name:
        subs = protein_desc.get("submissionNames", [])
        if subs:
            protein_name = subs[0].get("fullName", {}).get("value", "")

    gene_name = ""
    genes = entry.get("genes", [])
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")

    seq = entry.get("sequence", {})
    seq_len = seq.get("length", "")

    protein_existence = entry.get("proteinExistence", "")
    if isinstance(protein_existence, str) and ":" in protein_existence:
        protein_existence = protein_existence.split(":", 1)[0].strip()

    annotation_score = entry.get("annotationScore", "")

    return {
        "accession": accession,
        "protein_name": protein_name,
        "gene_name": gene_name,
        "organism": organism,
        "review_status": reviewed,
        "protein_existence": protein_existence,
        "annotation_score": annotation_score,
        "sequence_length": seq_len,
    }


def ncbi_search_assembly(organism: str) -> Optional[dict]:
    """
    NCBI Assembly search via Entrez esearch + esummary.
    """
    term = quote(f'"{organism}"[Organism] AND latest[filter]')
    search_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=assembly&term={term}&retmax=20&retmode=json"
    )
    search_data = http_get_json(search_url)
    ids = search_data.get("esearchresult", {}).get("idlist", [])
    if not ids:
        return None

    summary_url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
        f"?db=assembly&id={','.join(ids)}&retmode=json"
    )
    summary_data = http_get_json(summary_url)
    result = summary_data.get("result", {})
    uid_list = [uid for uid in result.get("uids", []) if uid != "uids"]
    if not uid_list:
        return None

    candidates = [result[uid] for uid in uid_list if uid in result]

    def score(record: dict) -> tuple[int, int, int]:
        status = record.get("assemblystatus", "") or ""
        refseq = bool(record.get("rsuid"))
        level_rank = {
            "Complete Genome": 4,
            "Chromosome": 3,
            "Scaffold": 2,
            "Contig": 1,
        }
        meta = record.get("meta")
        scaffold_n50 = 0
        if isinstance(meta, dict):
            try:
                scaffold_n50 = int(meta.get("scaffold_n50") or 0)
            except (TypeError, ValueError):
                scaffold_n50 = 0
        return (
            1 if refseq else 0,
            level_rank.get(status, 0),
            scaffold_n50,
        )

    candidates.sort(key=score, reverse=True)
    return candidates[0]


def parse_ncbi_assembly(summary: Optional[dict]) -> dict:
    if not summary:
        return {
            "assembly_accession": "",
            "assembly_level": "",
            "assembly_name": "",
            "refseq_or_genbank": "",
            "contig_n50": "",
            "scaffold_n50": "",
            "annotation_available": "",
            "assembly_status": "",
        }

    accession = summary.get("assemblyaccession", "")
    level = summary.get("assemblystatus", "")
    name = summary.get("assemblyname", "")
    refseq_or_genbank = "RefSeq" if summary.get("rsuid") else "GenBank"

    meta = summary.get("meta") if isinstance(summary.get("meta"), dict) else {}
    contig_n50 = meta.get("contig_n50", "")
    scaffold_n50 = meta.get("scaffold_n50", "")

    annotation_available = (
        "yes" if summary.get("gbuid") or summary.get("rsuid") else "unknown"
    )

    return {
        "assembly_accession": accession,
        "assembly_level": level,
        "assembly_name": name,
        "refseq_or_genbank": refseq_or_genbank,
        "contig_n50": contig_n50,
        "scaffold_n50": scaffold_n50,
        "annotation_available": annotation_available,
        "assembly_status": level,
    }


QUALITY_FLAG_STYLES = {
    "high": ("High", "#166534", "#dcfce7"),
    "medium": ("Medium", "#1e40af", "#dbeafe"),
    "review_protein_only": ("Protein only", "#92400e", "#fef3c7"),
    "low_or_manual_check": ("Low / check", "#991b1b", "#fee2e2"),
    "manual_check_no_accession": ("No accession", "#374151", "#f3f4f6"),
    "manual_check_uniprot_failed": ("UniProt failed", "#374151", "#f3f4f6"),
}


def write_html_report(rows: list[dict[str, str]], out_path: Path) -> None:
    flag_counts = Counter(row.get("quality_flag", "") for row in rows)
    family_counts = Counter(row.get("family", "") for row in rows)
    organisms = sorted({row.get("organism", "") for row in rows if row.get("organism")})

    def badge(flag: str) -> str:
        label, color, bg = QUALITY_FLAG_STYLES.get(
            flag, (flag.replace("_", " ").title(), "#374151", "#f3f4f6")
        )
        return (
            f'<span class="badge" style="color:{color};background:{bg};border:1px solid {color}22">'
            f"{html.escape(label)}</span>"
        )

    def cell(text: str, extra_class: str = "") -> str:
        value = html.escape(text or "—")
        cls = f' class="{extra_class}"' if extra_class else ""
        return f"<td{cls}>{value}</td>"

    table_rows = []
    for row in rows:
        accession = row.get("accession", "")
        accession_cell = (
            f'<a href="https://www.uniprot.org/uniprotkb/{html.escape(accession)}" '
            f'target="_blank" rel="noopener">{html.escape(accession)}</a>'
            if accession
            else "—"
        )
        assembly = row.get("assembly_accession", "")
        assembly_cell = (
            f'<a href="https://www.ncbi.nlm.nih.gov/assembly/{html.escape(assembly)}" '
            f'target="_blank" rel="noopener">{html.escape(assembly)}</a>'
            if assembly
            else '<span class="muted">not found</span>'
        )
        error = row.get("error", "").strip()
        error_html = (
            f'<div class="error-note">{html.escape(error)}</div>' if error else ""
        )
        flag = row.get("quality_flag", "")
        search_blob = " ".join(
            [
                row.get("organism", ""),
                row.get("gene_name", ""),
                row.get("protein_name", ""),
                accession,
                row.get("family", ""),
            ]
        ).lower()
        data_attrs = (
            f'data-quality="{html.escape(flag)}" '
            f'data-family="{html.escape(row.get("family", ""))}" '
            f'data-organism="{html.escape(row.get("organism", ""))}" '
            f'data-search="{html.escape(search_blob)}"'
        )
        table_rows.append(
            f"<tr {data_attrs}>"
            f"<td>{badge(flag)}</td>"
            f"{cell(row.get('family', ''), 'nowrap')}"
            f"{cell(row.get('organism', ''), 'organism')}"
            f"{cell(row.get('gene_name', ''), 'gene')}"
            f"{cell(row.get('protein_name', ''))}"
            f"<td class='mono'>{accession_cell}</td>"
            f"{cell(row.get('review_status', ''))}"
            f"{cell(str(row.get('annotation_score', '')), 'num')}"
            f"{cell(str(row.get('sequence_length', '')), 'num')}"
            f"{cell(row.get('assembly_level', ''), 'assembly')}"
            f"<td class='mono'>{assembly_cell}</td>"
            f"{cell(row.get('refseq_or_genbank', ''))}"
            f"{cell(row.get('assembly_name', ''), 'muted')}"
            f"<td>{error_html}</td>"
            "</tr>"
        )

    summary_cards = "".join(
        f'<div class="card"><div class="card-value">{count}</div>'
        f'<div class="card-label">{html.escape(label)}</div></div>'
        for label, count in sorted(flag_counts.items(), key=lambda x: -x[1])
    )

    family_options = "".join(
        f'<option value="{html.escape(f)}">{html.escape(f)} ({c})</option>'
        for f, c in sorted(family_counts.items())
    )

    quality_options = "".join(
        f'<option value="{html.escape(f)}">{html.escape(f)} ({c})</option>'
        for f, c in sorted(flag_counts.items())
    )

    document = f"""<!DOCTYPE html>
<html lang="en">
<head>
  <meta charset="utf-8">
  <meta name="viewport" content="width=device-width, initial-scale=1">
  <title>Gap junction quality overview</title>
  <style>
    :root {{
      --bg: #0f172a;
      --panel: #1e293b;
      --panel-2: #334155;
      --text: #f8fafc;
      --muted: #94a3b8;
      --accent: #38bdf8;
      --border: #475569;
    }}
    * {{ box-sizing: border-box; }}
    body {{
      margin: 0;
      font-family: "Segoe UI", system-ui, -apple-system, sans-serif;
      background: linear-gradient(160deg, #0f172a 0%, #1e293b 45%, #0f172a 100%);
      color: var(--text);
      min-height: 100vh;
    }}
    .wrap {{ max-width: 1600px; margin: 0 auto; padding: 2rem 1.5rem 3rem; }}
    h1 {{ margin: 0 0 0.35rem; font-size: 1.75rem; font-weight: 700; }}
    .subtitle {{ color: var(--muted); margin-bottom: 1.5rem; }}
    .stats {{
      display: grid;
      grid-template-columns: repeat(auto-fit, minmax(120px, 1fr));
      gap: 0.75rem;
      margin-bottom: 1.25rem;
    }}
    .card {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 0.9rem 1rem;
    }}
    .card-value {{ font-size: 1.6rem; font-weight: 700; color: var(--accent); }}
    .card-label {{ font-size: 0.8rem; color: var(--muted); margin-top: 0.2rem; word-break: break-word; }}
    .toolbar {{
      display: flex;
      flex-wrap: wrap;
      gap: 0.75rem;
      align-items: center;
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 12px;
      padding: 1rem;
      margin-bottom: 1rem;
    }}
    .toolbar label {{ font-size: 0.85rem; color: var(--muted); display: flex; flex-direction: column; gap: 0.35rem; }}
    input, select {{
      background: var(--panel-2);
      border: 1px solid var(--border);
      color: var(--text);
      border-radius: 8px;
      padding: 0.5rem 0.75rem;
      font-size: 0.95rem;
      min-width: 180px;
    }}
    input:focus, select:focus {{ outline: 2px solid var(--accent); outline-offset: 1px; }}
    .count-pill {{
      margin-left: auto;
      background: var(--panel-2);
      border-radius: 999px;
      padding: 0.45rem 0.9rem;
      font-size: 0.9rem;
      color: var(--muted);
    }}
    .table-shell {{
      background: var(--panel);
      border: 1px solid var(--border);
      border-radius: 12px;
      overflow: auto;
      box-shadow: 0 12px 40px rgba(0,0,0,0.25);
    }}
    table {{ width: 100%; border-collapse: collapse; font-size: 0.9rem; }}
    thead th {{
      position: sticky;
      top: 0;
      background: #111827;
      color: #e2e8f0;
      text-align: left;
      padding: 0.85rem 0.75rem;
      border-bottom: 2px solid var(--accent);
      white-space: nowrap;
      font-size: 0.78rem;
      text-transform: uppercase;
      letter-spacing: 0.04em;
    }}
    tbody td {{
      padding: 0.7rem 0.75rem;
      border-bottom: 1px solid var(--border);
      vertical-align: top;
    }}
    tbody tr:hover {{ background: rgba(56, 189, 248, 0.08); }}
    tbody tr.hidden {{ display: none; }}
    .badge {{
      display: inline-block;
      padding: 0.2rem 0.55rem;
      border-radius: 999px;
      font-size: 0.75rem;
      font-weight: 700;
      white-space: nowrap;
    }}
    .organism {{ font-weight: 600; min-width: 160px; }}
    .gene {{ font-weight: 600; color: #7dd3fc; }}
    .assembly {{ font-weight: 600; }}
    .mono {{ font-family: ui-monospace, monospace; font-size: 0.82rem; }}
    .num {{ text-align: center; }}
    .muted {{ color: var(--muted); font-size: 0.85rem; }}
    a {{ color: #7dd3fc; text-decoration: none; }}
    a:hover {{ text-decoration: underline; }}
    .error-note {{
      color: #fca5a5;
      font-size: 0.8rem;
      max-width: 220px;
    }}
    .legend {{
      display: flex;
      flex-wrap: wrap;
      gap: 0.5rem;
      margin: 0.75rem 0 1rem;
    }}
    .meta-line {{ color: var(--muted); font-size: 0.85rem; margin-bottom: 0.5rem; }}
  </style>
</head>
<body>
  <div class="wrap">
    <h1>Gap junction reference quality</h1>
    <p class="subtitle">
      {len(rows)} proteins · {len(organisms)} species ·
      generated from <code>gap_junction_quality_table.csv</code>
    </p>
    <div class="stats">
      <div class="card"><div class="card-value">{len(rows)}</div><div class="card-label">Proteins</div></div>
      <div class="card"><div class="card-value">{len(organisms)}</div><div class="card-label">Species</div></div>
      {summary_cards}
    </div>
    <div class="legend">
      <span class="badge" style="color:#166534;background:#dcfce7">High</span>
      <span class="badge" style="color:#1e40af;background:#dbeafe">Medium</span>
      <span class="badge" style="color:#92400e;background:#fef3c7">Protein only</span>
      <span class="badge" style="color:#991b1b;background:#fee2e2">Low / check</span>
    </div>
    <div class="toolbar">
      <label>Search
        <input id="search" type="search" placeholder="organism, gene, accession…">
      </label>
      <label>Family
        <select id="family-filter">
          <option value="">All families</option>
          {family_options}
        </select>
      </label>
      <label>Quality
        <select id="quality-filter">
          <option value="">All quality levels</option>
          {quality_options}
        </select>
      </label>
      <span class="count-pill" id="visible-count">{len(rows)} shown</span>
    </div>
    <p class="meta-line">Tip: click column headers to sort. Links open UniProt / NCBI Assembly.</p>
    <div class="table-shell">
      <table id="quality-table">
        <thead>
          <tr>
            <th data-sort="quality">Quality</th>
            <th data-sort="family">Family</th>
            <th data-sort="organism">Organism</th>
            <th data-sort="gene">Gene</th>
            <th data-sort="protein">Protein</th>
            <th data-sort="accession">UniProt</th>
            <th data-sort="review">Review</th>
            <th data-sort="score">Score</th>
            <th data-sort="length">Length</th>
            <th data-sort="assembly-level">Assembly</th>
            <th data-sort="assembly-acc">Accession</th>
            <th data-sort="source">Source</th>
            <th data-sort="assembly-name">Assembly name</th>
            <th>Notes</th>
          </tr>
        </thead>
        <tbody>
          {"".join(table_rows)}
        </tbody>
      </table>
    </div>
  </div>
  <script>
    const searchInput = document.getElementById("search");
    const familyFilter = document.getElementById("family-filter");
    const qualityFilter = document.getElementById("quality-filter");
    const visibleCount = document.getElementById("visible-count");
    const tbody = document.querySelector("#quality-table tbody");
    const rows = Array.from(tbody.querySelectorAll("tr"));

    function applyFilters() {{
      const q = searchInput.value.trim().toLowerCase();
      const family = familyFilter.value;
      const quality = qualityFilter.value;
      let shown = 0;
      rows.forEach((row) => {{
        const matchSearch = !q || row.dataset.search.includes(q);
        const matchFamily = !family || row.dataset.family === family;
        const matchQuality = !quality || row.dataset.quality === quality;
        const visible = matchSearch && matchFamily && matchQuality;
        row.classList.toggle("hidden", !visible);
        if (visible) shown += 1;
      }});
      visibleCount.textContent = shown + " shown";
    }}

    [searchInput, familyFilter, qualityFilter].forEach((el) =>
      el.addEventListener("input", applyFilters)
    );

    const colIndex = {{
      quality: 0, family: 1, organism: 2, gene: 3, protein: 4, accession: 5,
      review: 6, score: 7, length: 8, "assembly-level": 9, "assembly-acc": 10,
      source: 11, "assembly-name": 12
    }};
    let sortKey = "organism";
    let sortAsc = true;

    document.querySelectorAll("th[data-sort]").forEach((th) => {{
      th.style.cursor = "pointer";
      th.addEventListener("click", () => {{
        const key = th.dataset.sort;
        if (sortKey === key) sortAsc = !sortAsc;
        else {{ sortKey = key; sortAsc = true; }}
        const idx = colIndex[key];
        rows.sort((a, b) => {{
          const av = a.children[idx].innerText.trim().toLowerCase();
          const bv = b.children[idx].innerText.trim().toLowerCase();
          if (av < bv) return sortAsc ? -1 : 1;
          if (av > bv) return sortAsc ? 1 : -1;
          return 0;
        }});
        rows.forEach((row) => tbody.appendChild(row));
      }});
    }});
  </script>
</body>
</html>
"""
    out_path.write_text(document, encoding="utf-8")


def make_quality_flag(
    review_status: str, assembly_level: str, annotation_score: str
) -> str:
    level = (assembly_level or "").lower()

    try:
        ann_score = int(annotation_score) if annotation_score != "" else 0
    except (TypeError, ValueError):
        ann_score = 0

    if (
        review_status == "reviewed"
        and ("chromosome" in level or "complete" in level)
        and ann_score >= 4
    ):
        return "high"
    if review_status == "reviewed" and (
        "scaffold" in level or "chromosome" in level or "complete" in level
    ):
        return "medium"
    if review_status == "reviewed":
        return "review_protein_only"
    return "low_or_manual_check"


def main() -> None:
    fasta_files = find_fasta_files()
    if not fasta_files:
        print("No FASTA files found under", REFERENCES_DIR)
        sys.exit(1)

    print(f"Processing {len(fasta_files)} FASTA file(s)...")
    rows = []
    for index, (family, fasta_path) in enumerate(fasta_files, start=1):
        header = read_first_fasta_header(fasta_path)
        accession = extract_uniprot_accession(fasta_path, header)

        base_row = {
            "family": family,
            "file_path": str(fasta_path),
            "header": header,
            "accession": accession or "",
        }

        if not accession:
            rows.append(
                base_row
                | {
                    "protein_name": "",
                    "gene_name": "",
                    "organism": "",
                    "review_status": "",
                    "protein_existence": "",
                    "annotation_score": "",
                    "sequence_length": "",
                    "assembly_accession": "",
                    "assembly_level": "",
                    "assembly_name": "",
                    "refseq_or_genbank": "",
                    "contig_n50": "",
                    "scaffold_n50": "",
                    "annotation_available": "",
                    "assembly_status": "",
                    "quality_flag": "manual_check_no_accession",
                    "error": "Could not extract UniProt accession",
                }
            )
            continue

        print(f"[{index}/{len(fasta_files)}] {accession} ({fasta_path.name})")

        try:
            uni = parse_uniprot(get_uniprot_entry(accession))
            time.sleep(0.2)
        except Exception as exc:
            rows.append(
                base_row
                | {
                    "protein_name": "",
                    "gene_name": "",
                    "organism": "",
                    "review_status": "",
                    "protein_existence": "",
                    "annotation_score": "",
                    "sequence_length": "",
                    "assembly_accession": "",
                    "assembly_level": "",
                    "assembly_name": "",
                    "refseq_or_genbank": "",
                    "contig_n50": "",
                    "scaffold_n50": "",
                    "annotation_available": "",
                    "assembly_status": "",
                    "quality_flag": "manual_check_uniprot_failed",
                    "error": f"UniProt fetch failed: {exc}",
                }
            )
            continue

        error = ""
        try:
            ncbi_summary = ncbi_search_assembly(uni["organism"])
            ncbi = parse_ncbi_assembly(ncbi_summary)
            time.sleep(0.2)
        except Exception as exc:
            ncbi = parse_ncbi_assembly(None)
            error = f"NCBI assembly fetch failed: {exc}"

        quality_flag = make_quality_flag(
            review_status=uni["review_status"],
            assembly_level=ncbi["assembly_level"],
            annotation_score=str(uni["annotation_score"]),
        )

        rows.append(base_row | uni | ncbi | {"quality_flag": quality_flag, "error": error})

    out_csv = OUTPUT_DIR / "gap_junction_quality_table.csv"
    fieldnames = [
        "family",
        "file_path",
        "header",
        "accession",
        "protein_name",
        "gene_name",
        "organism",
        "review_status",
        "protein_existence",
        "annotation_score",
        "sequence_length",
        "assembly_accession",
        "assembly_level",
        "assembly_name",
        "refseq_or_genbank",
        "contig_n50",
        "scaffold_n50",
        "annotation_available",
        "assembly_status",
        "quality_flag",
        "error",
    ]

    with out_csv.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    flags = {}
    for row in rows:
        flags[row["quality_flag"]] = flags.get(row["quality_flag"], 0) + 1
    out_html = OUTPUT_DIR / "gap_junction_quality_table.html"
    write_html_report(rows, out_html)

    print(f"\nWrote: {out_csv}")
    print(f"Wrote: {out_html}")
    print("Quality flags:", flags)
    print("Open the HTML file in your browser for a visual overview.")


def html_only_from_csv(csv_path: Path) -> None:
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        print(f"No rows in {csv_path}")
        sys.exit(1)
    out_html = OUTPUT_DIR / "gap_junction_quality_table.html"
    write_html_report(rows, out_html)
    print(f"Wrote: {out_html} ({len(rows)} rows)")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build gap junction quality table.")
    parser.add_argument(
        "--html-only",
        action="store_true",
        help="Regenerate HTML report from existing CSV (no API calls).",
    )
    parser.add_argument(
        "--csv",
        type=Path,
        default=OUTPUT_DIR / "gap_junction_quality_table.csv",
        help="CSV path used with --html-only.",
    )
    return parser.parse_args()


if __name__ == "__main__":
    args = parse_args()
    if args.html_only:
        html_only_from_csv(args.csv)
    else:
        main()
