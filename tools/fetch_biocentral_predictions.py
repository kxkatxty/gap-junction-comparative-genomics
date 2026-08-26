from __future__ import annotations

import argparse
import hashlib
import json
import time
import urllib.error
import urllib.request
from pathlib import Path

import pandas as pd

from fetch_protein_features import (
    MAPPING_CSV,
    WHITELIST_CSV,
    load_mapping,
    load_whitelist,
    match_whitelist_entries,
)


PROJECT_ROOT = Path("project")
RESULTS_DIR = PROJECT_ROOT / "results" / "exon_structures"
BIOCENTRAL_CSV = RESULTS_DIR / "protein_structure_predictions.csv"
CACHE_DIR = RESULTS_DIR / "biocentral_cache"

BIOCENTRAL_BASE = "https://biocentral.rostlab.org"
PREDICT_URL = f"{BIOCENTRAL_BASE}/api/v1/prediction_service/predict"
TASK_STATUS_URL = f"{BIOCENTRAL_BASE}/api/v1/biocentral_service/task_status/{{task_id}}"

MODELS = ("ProtT5SecondaryStructure", "Seth", "TMbed")
USER_AGENT = "BachelorCursor/1.0 (academic research)"

POLL_INITIAL_WAIT_SEC = 60
POLL_INTERVAL_SEC = 15
POLL_MAX_WAIT_SEC = 900
BATCH_PAUSE_SEC = 25
BATCH_SIZE = 4

SETH_DISORDER_THRESHOLD = 8.0
MIN_REGION_LENGTH = 3

SEC_STRUCT_COLORS = {
    "helix": "#FF6B6B",
    "sheet": "#4DABF7",
    "coil": "#E9ECEF",
}
SEC_STRUCT_LABELS = {"H": "helix", "E": "sheet", "L": "coil"}
TM_POSITIVE = {
    "H": ("tm_helix", "TM helix", "#9B59B6"),
    "h": ("tm_helix", "TM helix", "#9B59B6"),
    "B": ("tm_beta", "TM beta", "#7D3C98"),
    "b": ("tm_beta", "TM beta", "#7D3C98"),
    "S": ("signal", "Signal", "#F1C40F"),
}


def read_fasta_sequence(path: Path) -> str:
    chunks: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                chunks.append(line.strip())
    return "".join(chunks)


def sequence_hash(sequence: str) -> str:
    return hashlib.sha256(sequence.encode("utf-8")).hexdigest()[:16]


def http_json(
    url: str,
    *,
    data: dict | None = None,
    method: str | None = None,
    timeout: int = 90,
) -> dict:
    headers = {"Accept": "application/json", "User-Agent": USER_AGENT}
    payload = None
    if data is not None:
        headers["Content-Type"] = "application/json"
        payload = json.dumps(data).encode("utf-8")
    request = urllib.request.Request(
        url,
        data=payload,
        headers=headers,
        method=method or ("POST" if payload is not None else "GET"),
    )
    with urllib.request.urlopen(request, timeout=timeout) as response:
        return json.loads(response.read().decode("utf-8"))


def submit_prediction_batch(sequence_input: dict[str, str]) -> str:
    payload = {"model_names": list(MODELS), "sequence_input": sequence_input}
    response = http_json(PREDICT_URL, data=payload, method="POST")
    return str(response["task_id"])


def poll_task(task_id: str) -> dict:
    url = TASK_STATUS_URL.format(task_id=task_id)
    time.sleep(POLL_INITIAL_WAIT_SEC)
    deadline = time.time() + POLL_MAX_WAIT_SEC
    while time.time() < deadline:
        try:
            response = http_json(url)
        except urllib.error.HTTPError as exc:
            if exc.code == 403:
                time.sleep(POLL_INTERVAL_SEC)
                continue
            raise
        for dto in response.get("dtos", []):
            status = dto.get("status")
            if status == "FINISHED":
                return dto
            if status == "FAILED":
                raise RuntimeError(f"Biocentral task failed: {dto}")
        time.sleep(POLL_INTERVAL_SEC)
    raise TimeoutError(f"Timed out waiting for Biocentral task {task_id}")


def merge_class_intervals(
    values: str,
    class_map: dict[str, tuple[str, str, str]],
    accession: str,
    track: str,
    source: str,
) -> list[dict]:
    features: list[dict] = []
    if not values:
        return features

    start: int | None = None
    current_key: str | None = None
    for index, char in enumerate(values, start=1):
        mapped = class_map.get(char)
        if mapped is None:
            if start is not None:
                feature_id, label, color = class_map[current_key]  # type: ignore[index]
                if index - start >= MIN_REGION_LENGTH:
                    features.append(
                        _feature_row(
                            accession,
                            track,
                            feature_id,
                            label,
                            start,
                            index - 1,
                            color,
                            source,
                        )
                    )
                start = None
                current_key = None
            continue

        feature_id, _, _ = mapped
        if start is None:
            start = index
            current_key = char
            continue

        prev_id, _, _ = class_map[current_key]  # type: ignore[index]
        if feature_id != prev_id:
            _, label, color = class_map[current_key]  # type: ignore[index]
            if index - start >= MIN_REGION_LENGTH:
                features.append(
                    _feature_row(
                        accession,
                        track,
                        prev_id,
                        label,
                        start,
                        index - 1,
                        color,
                        source,
                    )
                )
            start = index
            current_key = char

    if start is not None and current_key is not None:
        feature_id, label, color = class_map[current_key]
        end = len(values)
        if end - start + 1 >= MIN_REGION_LENGTH:
            features.append(
                _feature_row(
                    accession,
                    track,
                    feature_id,
                    label,
                    start,
                    end,
                    color,
                    source,
                )
            )
    return features


def merge_threshold_intervals(
    scores: list[float],
    accession: str,
    *,
    threshold: float,
    below_is_positive: bool = True,
) -> list[dict]:
    features: list[dict] = []
    if not scores:
        return features

    def is_positive(score: float) -> bool:
        return score < threshold if below_is_positive else score >= threshold

    start: int | None = None
    for index, score in enumerate(scores, start=1):
        active = is_positive(float(score))
        if active and start is None:
            start = index
        elif not active and start is not None:
            if index - start >= MIN_REGION_LENGTH:
                features.append(
                    _feature_row(
                        accession,
                        "disorder",
                        "disorder",
                        "Disordered",
                        start,
                        index - 1,
                        "#F39C12",
                        "biocentral_seth",
                    )
                )
            start = None

    if start is not None:
        end = len(scores)
        if end - start + 1 >= MIN_REGION_LENGTH:
            features.append(
                _feature_row(
                    accession,
                    "disorder",
                    "disorder",
                    "Disordered",
                    start,
                    end,
                    "#F39C12",
                    "biocentral_seth",
                )
            )
    return features


def _feature_row(
    accession: str,
    track: str,
    feature_id: str,
    feature_label: str,
    aa_start: int,
    aa_end: int,
    color: str,
    source: str,
) -> dict:
    return {
        "uniprot_accession": accession.upper(),
        "track": track,
        "feature_id": feature_id,
        "feature_label": feature_label,
        "aa_start": aa_start,
        "aa_end": aa_end,
        "color": color,
        "source": source,
    }


def parse_secondary_structure(accession: str, values: str) -> list[dict]:
    class_map = {
        key: (SEC_STRUCT_LABELS[key], SEC_STRUCT_LABELS[key].title(), SEC_STRUCT_COLORS[SEC_STRUCT_LABELS[key]])
        for key in SEC_STRUCT_LABELS
    }
    features = merge_class_intervals(
        values,
        class_map,
        accession,
        track="secondary_structure",
        source="biocentral_prott5",
    )
    return [feat for feat in features if feat["feature_id"] != "coil"]


def parse_transmembrane(accession: str, values: str) -> list[dict]:
    return merge_class_intervals(
        values,
        TM_POSITIVE,
        accession,
        track="transmembrane",
        source="biocentral_tmbed",
    )


def parse_seth_disorder(accession: str, value) -> list[dict]:
    if isinstance(value, list):
        scores = [float(v) for v in value]
    elif isinstance(value, str):
        scores = [float(part) for part in value.split(",") if part.strip()]
    else:
        return []
    return merge_threshold_intervals(
        scores,
        accession,
        threshold=SETH_DISORDER_THRESHOLD,
        below_is_positive=True,
    )


def predictions_to_features(accession: str, predictions: list[dict]) -> list[dict]:
    features: list[dict] = []
    for prediction in predictions:
        model_name = prediction.get("model_name", "")
        prediction_name = prediction.get("prediction_name", "")
        value = prediction.get("value")

        if model_name == "ProtT5SecondaryStructure" and prediction_name == "d3_Yhat":
            if isinstance(value, str):
                features.extend(parse_secondary_structure(accession, value))
        elif model_name == "TMbed" and prediction_name == "trans_membrane":
            if isinstance(value, str):
                features.extend(parse_transmembrane(accession, value))
        elif model_name == "Seth" and prediction_name == "disorder_chezod":
            features.extend(parse_seth_disorder(accession, value))
    return features


def load_targets(limit: int | None = None) -> pd.DataFrame:
    whitelist = load_whitelist(WHITELIST_CSV)
    mapping = load_mapping(MAPPING_CSV)
    targets = match_whitelist_entries(whitelist, mapping)
    if targets.empty:
        return targets

    rows: list[dict] = []
    for _, row in targets.iterrows():
        fasta_path = Path(str(row["fasta_path"]))
        if not fasta_path.exists():
            continue
        sequence = read_fasta_sequence(fasta_path)
        if not sequence:
            continue
        rows.append(
            {
                "plot_group": row["plot_group"],
                "organism": row["organism"],
                "uniprot_accession": str(row["uniprot_accession"]).upper(),
                "fasta_path": str(fasta_path),
                "sequence": sequence,
                "sequence_hash": sequence_hash(sequence),
            }
        )

    out = pd.DataFrame(rows).drop_duplicates(subset=["uniprot_accession", "sequence_hash"])
    out = out.sort_values(["plot_group", "organism", "uniprot_accession"])
    if limit is not None:
        out = out.head(limit)
    return out


def load_existing_predictions() -> tuple[pd.DataFrame, set[str]]:
    if not BIOCENTRAL_CSV.exists():
        return pd.DataFrame(), set()
    cached = pd.read_csv(BIOCENTRAL_CSV)
    if "sequence_hash" not in cached.columns:
        return cached, set()
    done = set(cached["sequence_hash"].astype(str).unique())
    return cached, done


def save_cache(accession: str, task_id: str, dto: dict, sequence_hash_value: str) -> None:
    CACHE_DIR.mkdir(parents=True, exist_ok=True)
    cache_file = CACHE_DIR / f"{accession}_{sequence_hash_value}.json"
    payload = {"task_id": task_id, "sequence_hash": sequence_hash_value, "dto": dto}
    cache_file.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def fetch_batch(batch: pd.DataFrame) -> list[dict]:
    sequence_input = {
        str(row["uniprot_accession"]): str(row["sequence"])
        for _, row in batch.iterrows()
    }
    task_id = submit_prediction_batch(sequence_input)
    accessions = ", ".join(sequence_input)
    print(f"  submitted {len(sequence_input)} sequence(s): {accessions}")
    print(f"  task_id={task_id}")

    dto = poll_task(task_id)
    rows: list[dict] = []
    for _, item in batch.iterrows():
        accession = str(item["uniprot_accession"])
        seq_hash = str(item["sequence_hash"])
        save_cache(accession, task_id, dto, seq_hash)
        predictions = dto.get("predictions", {}).get(accession, [])
        parsed = predictions_to_features(accession, predictions)
        for feature in parsed:
            feature["sequence_hash"] = seq_hash
        rows.extend(parsed)
        print(f"    {accession}: {len(parsed)} intervals")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Fetch Rostlab/Biocentral structure predictions for whitelist proteins."
    )
    parser.add_argument(
        "--limit",
        type=int,
        default=None,
        help="Fetch at most N unique sequences (useful for testing).",
    )
    parser.add_argument(
        "--accession",
        action="append",
        default=None,
        help="Fetch only these UniProt accessions (repeatable).",
    )
    parser.add_argument(
        "--force",
        action="store_true",
        help="Re-fetch even if sequence_hash is already cached.",
    )
    args = parser.parse_args()

    targets = load_targets(limit=args.limit)
    if args.accession:
        wanted = {acc.strip().upper() for acc in args.accession}
        targets = targets[targets["uniprot_accession"].astype(str).str.upper().isin(wanted)]
    if targets.empty:
        print("No whitelist proteins with FASTA sequences found.")
        return

    cached_df, cached_hashes = load_existing_predictions()
    if args.force:
        pending = targets
        base_rows: list[dict] = []
    else:
        pending = targets[~targets["sequence_hash"].astype(str).isin(cached_hashes)]
        base_rows = (
            cached_df.to_dict("records") if not cached_df.empty else []
        )

    print(
        f"Targets: {len(targets)} unique sequences; "
        f"pending fetch: {len(pending)}; cached: {len(cached_hashes)}."
    )
    if pending.empty:
        print(f"Nothing to fetch. Existing predictions: {BIOCENTRAL_CSV}")
        return

    fetched_rows: list[dict] = []
    batches = [
        pending.iloc[start : start + BATCH_SIZE]
        for start in range(0, len(pending), BATCH_SIZE)
    ]
    for batch_index, batch in enumerate(batches, start=1):
        print(f"Batch {batch_index}/{len(batches)}")
        try:
            fetched_rows.extend(fetch_batch(batch))
        except Exception as exc:  # noqa: BLE001 - keep pipeline running per batch
            print(f"  warning: batch failed: {exc}")
        if batch_index < len(batches):
            time.sleep(BATCH_PAUSE_SEC)

    if not fetched_rows and not base_rows:
        print("No predictions fetched.")
        return

    out = pd.DataFrame(base_rows + fetched_rows)
    dedupe_cols = [
        "uniprot_accession",
        "sequence_hash",
        "track",
        "feature_id",
        "aa_start",
        "aa_end",
    ]
    dedupe_cols = [col for col in dedupe_cols if col in out.columns]
    out = out.drop_duplicates(subset=dedupe_cols)
    sort_cols = [col for col in ["uniprot_accession", "track", "aa_start"] if col in out.columns]
    out = out.sort_values(sort_cols)
    BIOCENTRAL_CSV.parent.mkdir(parents=True, exist_ok=True)
    out.to_csv(BIOCENTRAL_CSV, index=False)
    print(f"Saved {len(out)} structure intervals to {BIOCENTRAL_CSV}")


if __name__ == "__main__":
    main()
