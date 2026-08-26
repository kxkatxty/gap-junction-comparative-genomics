from __future__ import annotations

import csv
import json
import re
import sys
import time
from pathlib import Path
from urllib.parse import quote
from urllib.request import Request, urlopen
from urllib.error import HTTPError, URLError


PROJECT_ROOT = Path("project")
REFERENCES_DIR = PROJECT_ROOT / "data" / "references"
METADATA_DIR = PROJECT_ROOT / "metadata"
USER_AGENT = "gap-junction-search-downloader/1.0"

# UniProt may still index proteins under older species names.
UNIPROT_ORGANISM_ALIASES: dict[str, list[str]] = {
    "Chlorocebus sabaeus": ["Chlorocebus sabaeus", "Chlorocebus aethiops"],
    "Hylobates moloch": ["Hylobates moloch", "Hylobates lar"],
}


def ensure_directories() -> None:
    for path in [
        REFERENCES_DIR / "innexins",
        REFERENCES_DIR / "connexins",
        REFERENCES_DIR / "pannexins",
        METADATA_DIR,
    ]:
        path.mkdir(parents=True, exist_ok=True)


def sanitize_filename(text: str) -> str:
    text = text.strip().replace(" ", "_")
    return re.sub(r"[^A-Za-z0-9._-]", "_", text)


def family_to_folder(family: str) -> Path:
    family_norm = family.lower().strip()
    if family_norm in {"innexin", "innexins"}:
        return REFERENCES_DIR / "innexins"
    if family_norm in {"connexin", "connexins"}:
        return REFERENCES_DIR / "connexins"
    if family_norm in {"pannexin", "pannexins"}:
        return REFERENCES_DIR / "pannexins"
    raise ValueError(f"Unsupported family: {family}")


def search_uniprot(query: str, fields: str | None = None, size: int = 500) -> dict:
    """
    Search UniProtKB and return JSON.
    """
    base_url = "https://rest.uniprot.org/uniprotkb/search"
    params = [f"query={quote(query)}", "format=json", f"size={size}"]
    if fields:
        params.append(f"fields={quote(fields)}")
    url = f"{base_url}?{'&'.join(params)}"

    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_fasta(accession: str) -> str:
    url = f"https://rest.uniprot.org/uniprotkb/{accession}.fasta"
    request = Request(url, headers={"User-Agent": USER_AGENT})
    with urlopen(request, timeout=30) as response:
        text = response.read().decode("utf-8")
    if not text.startswith(">"):
        raise ValueError(f"{accession} did not return FASTA content.")
    return text


def uniprot_organism_names(organism: str) -> list[str]:
    names = UNIPROT_ORGANISM_ALIASES.get(organism, [organism])
    seen: set[str] = set()
    unique: list[str] = []
    for name in names:
        key = name.casefold()
        if key not in seen:
            seen.add(key)
            unique.append(name)
    return unique


def build_query(family: str, organism: str, reviewed_only: bool = True) -> str:
    """
    Build a UniProt query.
    """
    terms = [
        f"protein_name:{family}",
        f'organism_name:"{organism}"',
    ]
    if reviewed_only:
        terms.append("reviewed:true")
    return " AND ".join(terms)


def extract_entry_info(entry: dict) -> dict[str, str]:
    accession = entry.get("primaryAccession", "")
    protein_desc = entry.get("proteinDescription", {})
    recommended = protein_desc.get("recommendedName", {})
    full_name = recommended.get("fullName", {}).get("value", "")

    gene_name = ""
    genes = entry.get("genes", [])
    if genes:
        gene_name = genes[0].get("geneName", {}).get("value", "")

    organism = entry.get("organism", {}).get("scientificName", "")
    reviewed = (
        "reviewed"
        if entry.get("entryType", "").startswith("UniProtKB reviewed")
        else "unreviewed"
    )

    return {
        "accession": accession,
        "gene_name": gene_name,
        "protein_name": full_name,
        "organism": organism,
        "review_status": reviewed,
    }


def save_fasta(
    folder: Path, organism: str, gene_name: str, accession: str, fasta_text: str
) -> Path:
    organism_safe = sanitize_filename(organism)
    gene_safe = sanitize_filename(gene_name or accession)
    out_dir = folder / organism_safe
    out_dir.mkdir(parents=True, exist_ok=True)
    out_path = out_dir / f"{gene_safe}__{accession}.fasta"
    out_path.write_text(fasta_text, encoding="utf-8")
    return out_path


def main(config_csv: str) -> None:
    """
    config_csv columns:
    family,organism,reviewed_only

    Example:
    innexin,Drosophila melanogaster,true
    innexin,Caenorhabditis elegans,true
    connexin,Homo sapiens,true
    """
    ensure_directories()

    rows = []
    with Path(config_csv).open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"family", "organism", "reviewed_only"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing columns in config CSV: {', '.join(sorted(missing))}")
        for row in reader:
            rows.append({k: (v or "").strip() for k, v in row.items()})

    report = []

    for row in rows:
        family = row["family"]
        organism = row["organism"]
        reviewed_only = row["reviewed_only"].lower() == "true"

        folder = family_to_folder(family)
        entries: list[dict] = []
        used_query = ""

        try:
            for search_organism in uniprot_organism_names(organism):
                query = build_query(
                    family, search_organism, reviewed_only=reviewed_only
                )
                print(f"\nSearching UniProt for: {query}")
                results = search_uniprot(query, size=500)
                entries = results.get("results", [])
                if entries:
                    used_query = query
                    if search_organism != organism:
                        print(f"  Found via UniProt alias: {search_organism}")
                    break
                time.sleep(0.2)

            if not entries:
                print(f"  No results found for {family} in {organism}")
                report.append(
                    {
                        "family": family,
                        "organism": organism,
                        "accession": "",
                        "gene_name": "",
                        "protein_name": "",
                        "review_status": "",
                        "status": "no_results",
                        "file_path": "",
                        "error": "",
                    }
                )
                continue

            for entry in entries:
                info = extract_entry_info(entry)
                accession = info["accession"]

                try:
                    fasta_text = fetch_fasta(accession)
                    saved_path = save_fasta(
                        folder=folder,
                        organism=organism,
                        gene_name=info["gene_name"],
                        accession=accession,
                        fasta_text=fasta_text,
                    )

                    print(f"  [OK] {accession} {info['gene_name']} -> {saved_path}")
                    report.append(
                        {
                            "family": family,
                            "organism": organism,
                            "accession": accession,
                            "gene_name": info["gene_name"],
                            "protein_name": info["protein_name"],
                            "review_status": info["review_status"],
                            "status": "downloaded",
                            "file_path": str(saved_path),
                            "error": "",
                        }
                    )

                    time.sleep(0.2)

                except Exception as exc:
                    print(f"  [FAIL] {accession}: {exc}")
                    report.append(
                        {
                            "family": family,
                            "organism": organism,
                            "accession": accession,
                            "gene_name": info["gene_name"],
                            "protein_name": info["protein_name"],
                            "review_status": info["review_status"],
                            "status": "failed",
                            "file_path": "",
                            "error": str(exc),
                        }
                    )

        except (HTTPError, URLError, Exception) as exc:
            print(f"Search failed for {family} / {organism}: {exc}")
            used_query = used_query or f"{family} / {organism}"
            report.append(
                {
                    "family": family,
                    "organism": organism,
                    "accession": "",
                    "gene_name": "",
                    "protein_name": "",
                    "review_status": "",
                    "status": "search_failed",
                    "file_path": "",
                    "error": str(exc),
                }
            )

    report_path = METADATA_DIR / "uniprot_family_download_report.csv"
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        fieldnames = [
            "family",
            "organism",
            "accession",
            "gene_name",
            "protein_name",
            "review_status",
            "status",
            "file_path",
            "error",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report)

    print(f"\nReport written to: {report_path}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        print("Usage: python download_by_family.py <config.csv>")
        sys.exit(1)
    main(sys.argv[1])
