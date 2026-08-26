from __future__ import annotations

import argparse
import csv
import gzip
import shutil
import sys
import time
import xml.etree.ElementTree as ET
from pathlib import Path
from urllib.error import HTTPError, URLError
from urllib.parse import quote
from urllib.request import Request, urlopen


USER_AGENT = "gap-junction-gff3-downloader/1.0"

PROJECT_ROOT = Path("project")
ANNOTATION_DIR = PROJECT_ROOT / "data" / "annotations"
METADATA_DIR = PROJECT_ROOT / "metadata"

ANNOTATION_DIR.mkdir(parents=True, exist_ok=True)
METADATA_DIR.mkdir(parents=True, exist_ok=True)


def sanitize_filename(text: str) -> str:
    return (
        text.strip()
        .replace(" ", "_")
        .replace("/", "_")
        .replace("(", "")
        .replace(")", "")
    )


def read_species_csv(csv_path: Path) -> list[dict[str, str]]:
    with csv_path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle)
        required = {"family", "organism", "reviewed_only"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Missing required columns: {', '.join(sorted(missing))}")
        return [{k: (v or "").strip() for k, v in row.items()} for row in reader]


def http_get_text(url: str, timeout_s: int = 60, retries: int = 3) -> str:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            req = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(req, timeout=timeout_s) as resp:
                return resp.read().decode("utf-8")
        except HTTPError as exc:
            last_error = exc
            if exc.code == 429 and attempt < retries - 1:
                time.sleep(2.0 * (attempt + 1))
                continue
            raise
        except URLError as exc:
            last_error = exc
            if attempt < retries - 1:
                time.sleep(1.0 * (attempt + 1))
                continue
            raise
    raise last_error or RuntimeError("http_get_text failed")


def url_exists(url: str) -> bool:
    req = Request(url, method="HEAD", headers={"User-Agent": USER_AGENT})
    try:
        with urlopen(req, timeout=30) as resp:
            return 200 <= resp.status < 300
    except HTTPError as exc:
        if exc.code == 404:
            return False
        raise
    except URLError:
        return False


def entrez_esearch_assembly(organism: str, term: str, retmax: int = 25) -> list[str]:
    url = (
        "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
        f"?db=assembly&term={quote(term)}&retmode=xml&retmax={retmax}"
    )
    xml_text = http_get_text(url)
    root = ET.fromstring(xml_text)
    return [elem.text for elem in root.findall(".//IdList/Id") if elem.text]


def entrez_esummary_assembly(uids: list[str]) -> list[dict[str, str]]:
    if not uids:
        return []

    records: list[dict[str, str]] = []
    chunk_size = 20
    for start in range(0, len(uids), chunk_size):
        chunk = uids[start : start + chunk_size]
        url = (
            "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi"
            f"?db=assembly&id={','.join(chunk)}&retmode=xml"
        )
        xml_text = http_get_text(url)
        root = ET.fromstring(xml_text)
        for doc in root.findall(".//DocumentSummary"):
            rec: dict[str, str] = {}
            for child in list(doc):
                rec[child.tag] = (child.text or "").strip()
            records.append(rec)
        time.sleep(0.34)
    return records


def assembly_search_terms(organism: str) -> list[str]:
    quoted = f'"{organism}"[Organism]'
    return [
        f"{quoted} AND latest_refseq[filter]",
        f"{quoted} AND reference genome[Properties]",
        f"{quoted} AND representative genome[Properties]",
        f"{quoted} AND latest[filter]",
        f"{quoted}",
    ]


def fetch_assembly_candidates(organism: str) -> list[dict[str, str]]:
    seen_uids: set[str] = set()
    candidates: list[dict[str, str]] = []

    for term in assembly_search_terms(organism):
        uids = entrez_esearch_assembly(organism, term)
        new_uids = [uid for uid in uids if uid not in seen_uids]
        if not new_uids:
            continue
        seen_uids.update(new_uids)
        candidates.extend(entrez_esummary_assembly(new_uids))
        time.sleep(0.34)

    return candidates


def rank_assemblies(records: list[dict[str, str]]) -> list[dict[str, str]]:
    def score(rec: dict[str, str]) -> tuple[int, int, int, int]:
        ftp_refseq = bool(rec.get("FtpPath_RefSeq"))
        category = (rec.get("RefSeq_category", "") or "").lower()
        level = (rec.get("AssemblyStatus", "") or "").lower()
        accession = rec.get("AssemblyAccession", "")

        refseq_score = 2 if ftp_refseq else 0
        if accession.startswith("GCF_"):
            refseq_score += 1

        if "reference genome" in category:
            category_score = 3
        elif "representative genome" in category:
            category_score = 2
        elif category:
            category_score = 1
        else:
            category_score = 0

        if "complete genome" in level:
            level_score = 4
        elif "chromosome" in level:
            level_score = 3
        elif "scaffold" in level:
            level_score = 2
        elif "contig" in level:
            level_score = 1
        else:
            level_score = 0

        # Prefer annotated RefSeq reference builds over raw GenBank-only assemblies.
        genbank_only = 0 if ftp_refseq else -2

        return (refseq_score, category_score, level_score, genbank_only)

    unique: dict[str, dict[str, str]] = {}
    for rec in records:
        accession = rec.get("AssemblyAccession", "")
        if not accession:
            continue
        if accession not in unique or score(rec) > score(unique[accession]):
            unique[accession] = rec

    return sorted(unique.values(), key=score, reverse=True)


def gff_candidate_urls(assembly_record: dict[str, str]) -> list[tuple[str, str, str]]:
    """
    Return list of (url, assembly_base, source) to try.
    Prefer RefSeq FTP path because GFF is usually published there.
    """
    options: list[tuple[str, str, str]] = []
    for source_key, label in (
        ("FtpPath_RefSeq", "RefSeq"),
        ("FtpPath_GenBank", "GenBank"),
    ):
        ftp_path = assembly_record.get(source_key, "")
        if not ftp_path:
            continue
        assembly_base = ftp_path.rstrip("/").split("/")[-1]
        https_path = ftp_path.replace("ftp://", "https://")
        options.append(
            (f"{https_path}/{assembly_base}_genomic.gff.gz", assembly_base, label)
        )
        options.append(
            (f"{https_path}/{assembly_base}_genomic.gff", assembly_base, label)
        )
    return options


def existing_gff_files(organism_dir: Path) -> list[Path]:
    return sorted(
        list(organism_dir.glob("*.gff"))
        + list(organism_dir.glob("*.gff.gz"))
        + list(organism_dir.glob("*.gff3"))
        + list(organism_dir.glob("*.gff3.gz"))
    )


def _validate_downloaded_file(path: Path) -> None:
    if not path.exists() or path.stat().st_size < 1000:
        raise OSError(f"download too small: {path}")
    with path.open("rb") as handle:
        if handle.read(2) == b"\x1f\x8b":
            with gzip.open(path, "rb") as gz_handle:
                while gz_handle.read(1024 * 1024):
                    pass


def download_file(url: str, out_path: Path, retries: int = 5, timeout_s: int = 900) -> None:
    last_error: Exception | None = None
    for attempt in range(retries):
        try:
            out_path.parent.mkdir(parents=True, exist_ok=True)
            req = Request(url, headers={"User-Agent": USER_AGENT})
            with urlopen(req, timeout=timeout_s) as resp, out_path.open("wb") as out_f:
                shutil.copyfileobj(resp, out_f)
            _validate_downloaded_file(out_path)
            return
        except HTTPError as exc:
            last_error = exc
            out_path.unlink(missing_ok=True)
            if exc.code == 429 and attempt < retries - 1:
                time.sleep(5.0 * (attempt + 1))
                continue
            raise
        except (URLError, TimeoutError, OSError, gzip.BadGzipFile) as exc:
            last_error = exc
            out_path.unlink(missing_ok=True)
            if attempt < retries - 1:
                time.sleep(5.0 * (attempt + 1))
                continue
            raise
    raise last_error or RuntimeError("download_file failed")


def maybe_decompress_gzip(gz_path: Path, keep_gz: bool = True) -> Path:
    out_path = gz_path.with_suffix("")
    with gzip.open(gz_path, "rb") as f_in, out_path.open("wb") as f_out:
        shutil.copyfileobj(f_in, f_out)
    if not keep_gz:
        gz_path.unlink(missing_ok=True)
    return out_path


def download_gff_for_organism(
    organism: str,
    organism_dir: Path,
    decompress: bool,
) -> dict[str, str]:
    existing = existing_gff_files(organism_dir)
    if existing:
        path = existing[0]
        return {
            "assembly_name": "",
            "assembly_accession": "",
            "assembly_level": "",
            "refseq_category": "",
            "ftp_source": "",
            "gff3_url": "",
            "downloaded_gff_gz": str(path) if path.suffix == ".gz" else "",
            "downloaded_gff3": str(path.with_suffix("")) if path.suffix == ".gz" else str(path),
            "status": "skipped_existing",
            "error": "",
        }

    candidates = rank_assemblies(fetch_assembly_candidates(organism))
    if not candidates:
        raise ValueError("No assembly candidates found.")

    errors: list[str] = []
    for assembly in candidates[:12]:
        for gff_url, assembly_base, source in gff_candidate_urls(assembly):
            if not url_exists(gff_url):
                continue

            if gff_url.endswith(".gz"):
                out_path = organism_dir / f"{assembly_base}_genomic.gff.gz"
            else:
                out_path = organism_dir / f"{assembly_base}_genomic.gff"

            print(f"  Trying {assembly.get('AssemblyAccession')} ({source})")
            print(f"  Downloading: {gff_url}")
            download_file(gff_url, out_path)

            result = {
                "assembly_name": assembly.get("AssemblyName", ""),
                "assembly_accession": assembly.get("AssemblyAccession", ""),
                "assembly_level": assembly.get("AssemblyStatus", ""),
                "refseq_category": assembly.get("RefSeq_category", ""),
                "ftp_source": assembly.get("FtpPath_RefSeq")
                or assembly.get("FtpPath_GenBank")
                or "",
                "gff3_url": gff_url,
                "downloaded_gff_gz": "",
                "downloaded_gff3": "",
                "status": "downloaded",
                "error": "",
            }

            if out_path.suffix == ".gz":
                result["downloaded_gff_gz"] = str(out_path)
                if decompress:
                    gff_path = maybe_decompress_gzip(out_path, keep_gz=True)
                    result["downloaded_gff3"] = str(gff_path)
            else:
                result["downloaded_gff3"] = str(out_path)

            return result

        errors.append(
            f"{assembly.get('AssemblyAccession', '?')}: no downloadable GFF URL"
        )

    raise ValueError("; ".join(errors[:5]))


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Download NCBI assembly genomic GFF for a species list."
    )
    parser.add_argument("csv", help="CSV with columns family,organism,reviewed_only")
    parser.add_argument(
        "--no-decompress",
        action="store_true",
        help="Keep only .gff.gz (do not unzip).",
    )
    parser.add_argument(
        "--sleep-s",
        type=float,
        default=0.5,
        help="Sleep between organisms to reduce NCBI rate limits.",
    )
    return parser.parse_args(argv)


def main(csv_file: str, decompress: bool, sleep_s: float) -> None:
    rows = read_species_csv(Path(csv_file))
    report_rows: list[dict[str, str]] = []

    for index, row in enumerate(rows, start=1):
        family = row["family"]
        organism = row["organism"]
        reviewed_only = row["reviewed_only"]

        family_dir = ANNOTATION_DIR / sanitize_filename(family.lower())
        family_dir.mkdir(parents=True, exist_ok=True)
        organism_dir = family_dir / sanitize_filename(organism)
        organism_dir.mkdir(parents=True, exist_ok=True)

        result: dict[str, str] = {
            "family": family,
            "organism": organism,
            "reviewed_only": reviewed_only,
            "assembly_name": "",
            "assembly_accession": "",
            "assembly_level": "",
            "refseq_category": "",
            "ftp_source": "",
            "gff3_url": "",
            "downloaded_gff_gz": "",
            "downloaded_gff3": "",
            "status": "",
            "error": "",
        }

        print(f"\n[{index}/{len(rows)}] {organism}")
        try:
            download_result = download_gff_for_organism(
                organism=organism,
                organism_dir=organism_dir,
                decompress=decompress,
            )
            result.update(download_result)
            if result["status"] == "skipped_existing":
                print(f"[SKIP] already present -> {result['downloaded_gff3'] or result['downloaded_gff_gz']}")
            else:
                print(
                    f"[OK] {organism} -> "
                    f"{result['downloaded_gff3'] or result['downloaded_gff_gz']}"
                )
        except (HTTPError, URLError, Exception) as exc:
            result["status"] = "failed"
            result["error"] = str(exc)
            print(f"[FAIL] {organism}: {exc}")

        report_rows.append(result)
        time.sleep(max(sleep_s, 0.0))

    report_path = METADATA_DIR / "gff3_download_report.csv"
    fieldnames = [
        "family",
        "organism",
        "reviewed_only",
        "assembly_name",
        "assembly_accession",
        "assembly_level",
        "refseq_category",
        "ftp_source",
        "gff3_url",
        "downloaded_gff_gz",
        "downloaded_gff3",
        "status",
        "error",
    ]
    with report_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(report_rows)

    ok = sum(1 for r in report_rows if r["status"] in {"downloaded", "skipped_existing"})
    failed = sum(1 for r in report_rows if r["status"] == "failed")
    print(f"\nReport written to: {report_path}")
    print(f"Success: {ok}/{len(report_rows)} | Failed: {failed}")


if __name__ == "__main__":
    args = parse_args()
    main(args.csv, decompress=not args.no_decompress, sleep_s=args.sleep_s)
