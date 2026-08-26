from __future__ import annotations

import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd


PROJECT_ROOT = Path("project")
METADATA_DIR = PROJECT_ROOT / "metadata"
RESULTS_DIR = PROJECT_ROOT / "results" / "exon_structures"
WHITELIST_CSV = METADATA_DIR / "gene_plot_whitelist.csv"
MAPPING_CSV = RESULTS_DIR / "protein_mapping.csv"
FEATURES_CSV = RESULTS_DIR / "protein_features.csv"

INTERPRO_PFAM_URL = "https://www.ebi.ac.uk/interpro/api/entry/pfam/protein/uniprot/{accession}"
UNIPROT_URL = "https://rest.uniprot.org/uniprotkb/{accession}.json"
REQUEST_DELAY_SEC = 0.25

PFAM_COLORS = ("#E67E22", "#8E44AD", "#2980B9", "#C0392B")


def load_whitelist(path: Path) -> pd.DataFrame:
    df = pd.read_csv(path)
    required = {"plot_group", "family", "reference_gene_symbol", "enabled"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"{path} missing columns: {sorted(missing)}")
    return df


def load_mapping(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(f"Missing {path}")
    return pd.read_csv(path)


def _is_enabled(value) -> bool:
    if pd.isna(value):
        return False
    return str(value).strip().lower() in {"1", "true", "yes", "y"}


def match_whitelist_entries(whitelist: pd.DataFrame, mapping: pd.DataFrame) -> pd.DataFrame:
    rows: list[dict] = []
    active = whitelist[whitelist["enabled"].map(_is_enabled)]

    for _, wl in active.iterrows():
        mask = mapping["family"].astype(str).str.lower() == str(wl["family"]).lower()
        symbol = str(wl["reference_gene_symbol"]).strip()
        mask &= mapping["reference_gene_symbol"].astype(str).str.lower() == symbol.lower()

        accession = wl.get("uniprot_accession", "")
        if pd.notna(accession) and str(accession).strip():
            mask &= (
                mapping["uniprot_accession"].astype(str).str.upper()
                == str(accession).strip().upper()
            )

        organism = wl.get("organism", "")
        if pd.notna(organism) and str(organism).strip():
            mask &= mapping["organism"].astype(str) == str(organism).strip()

        hits = mapping[mask]
        for _, hit in hits.iterrows():
            rows.append(
                {
                    "plot_group": wl["plot_group"],
                    "title_label": wl.get("title_label") or symbol,
                    "family": hit["family"],
                    "organism": hit["organism"],
                    "fasta_path": hit["fasta_path"],
                    "uniprot_accession": hit["uniprot_accession"],
                    "reference_gene_symbol": hit["reference_gene_symbol"],
                }
            )

    if not rows:
        return pd.DataFrame()

    out = pd.DataFrame(rows).drop_duplicates(
        subset=["plot_group", "fasta_path", "uniprot_accession"]
    )
    return out.sort_values(["plot_group", "organism", "reference_gene_symbol"])


def http_get_json(url: str) -> dict:
    request = urllib.request.Request(
        url,
        headers={"Accept": "application/json", "User-Agent": "BachelorCursor/1.0"},
    )
    with urllib.request.urlopen(request, timeout=60) as response:
        return json.loads(response.read().decode("utf-8"))


def fetch_pfam_features(accession: str) -> list[dict]:
    url = INTERPRO_PFAM_URL.format(accession=accession)
    try:
        payload = http_get_json(url)
    except urllib.error.HTTPError:
        return []

    features: list[dict] = []
    color_index = 0
    for entry in payload.get("results", []):
        meta = entry.get("metadata", {})
        label = meta.get("name") or meta.get("accession", "Pfam")
        pfam_id = meta.get("accession", "")
        display = f"{label}" if not pfam_id else f"{label}"
        color = PFAM_COLORS[color_index % len(PFAM_COLORS)]
        color_index += 1

        for protein in entry.get("proteins", []):
            for location in protein.get("entry_protein_locations", []):
                for fragment in location.get("fragments", []):
                    features.append(
                        {
                            "uniprot_accession": accession.upper(),
                            "track": "pfam",
                            "feature_id": pfam_id,
                            "feature_label": display,
                            "aa_start": int(fragment["start"]),
                            "aa_end": int(fragment["end"]),
                            "color": color,
                            "source": "interpro",
                        }
                    )
    return features


def fetch_uniprot_features(accession: str) -> list[dict]:
    url = UNIPROT_URL.format(accession=accession)
    try:
        payload = http_get_json(url)
    except urllib.error.HTTPError:
        return []

    features: list[dict] = []
    for feature in payload.get("features", []):
        ftype = feature.get("type", "")
        description = (feature.get("description") or "").strip()
        start = int(feature["location"]["start"]["value"])
        end = int(feature["location"]["end"]["value"])

        if ftype == "Coiled coil":
            features.append(
                {
                    "uniprot_accession": accession.upper(),
                    "track": "coiled_coil",
                    "feature_id": "coiled_coil",
                    "feature_label": "Coiled coil",
                    "aa_start": start,
                    "aa_end": end,
                    "color": "#27AE60",
                    "source": "uniprot",
                }
            )
        elif ftype == "Region" and "disorder" in description.lower():
            features.append(
                {
                    "uniprot_accession": accession.upper(),
                    "track": "disorder",
                    "feature_id": "disorder",
                    "feature_label": "Disordered",
                    "aa_start": start,
                    "aa_end": end,
                    "color": "#BDC3C7",
                    "source": "uniprot",
                }
            )
    return features


def fetch_features_for_accession(accession: str) -> list[dict]:
    accession = accession.upper()
    features = fetch_pfam_features(accession)
    time.sleep(REQUEST_DELAY_SEC)
    features.extend(fetch_uniprot_features(accession))
    time.sleep(REQUEST_DELAY_SEC)
    return features


def main() -> None:
    whitelist = load_whitelist(WHITELIST_CSV)
    mapping = load_mapping(MAPPING_CSV)
    targets = match_whitelist_entries(whitelist, mapping)

    if targets.empty:
        print("No proteins matched the whitelist.")
        return

    existing: set[str] = set()
    rows: list[dict] = []
    if FEATURES_CSV.exists():
        cached = pd.read_csv(FEATURES_CSV)
        rows.extend(cached.to_dict("records"))
        existing = set(cached["uniprot_accession"].astype(str).str.upper())

    accessions = sorted(targets["uniprot_accession"].astype(str).str.upper().unique())
    to_fetch = [acc for acc in accessions if acc not in existing]

    print(f"Whitelist matched {len(targets)} proteins ({len(accessions)} accessions).")
    print(f"Fetching {len(to_fetch)} new accessions...")

    for accession in to_fetch:
        print(f"  {accession}")
        try:
            rows.extend(fetch_features_for_accession(accession))
        except Exception as exc:  # noqa: BLE001 - keep pipeline running per accession
            print(f"    warning: failed for {accession}: {exc}")

    if not rows:
        print("No features fetched.")
        return

    out = pd.DataFrame(rows).drop_duplicates(
        subset=["uniprot_accession", "track", "feature_id", "aa_start", "aa_end"]
    )
    out = out.sort_values(["uniprot_accession", "track", "aa_start"])
    FEATURES_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(FEATURES_CSV, index=False)
    print(f"Saved {len(out)} feature intervals to {FEATURES_CSV}")


if __name__ == "__main__":
    main()
