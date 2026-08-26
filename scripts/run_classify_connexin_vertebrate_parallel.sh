#!/usr/bin/env bash
# Classify connexin annotation for vertebrate clades in parallel (5 jobs).
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
BATCH_DIR="${ROOT}/project/metadata/connexin_vertebrate_batches"
OUT_DIR="${ROOT}/project/metadata/connexin_vertebrate_classify_parts"
LOG_DIR="${ROOT}/project/results/connexin_vertebrate_classify/logs"
MERGED_STATUS="${ROOT}/project/metadata/family_annotation_status_connexin_vertebrate.csv"
MERGED_SUMMARY="${ROOT}/project/metadata/family_annotation_summary_connexin_vertebrate.csv"
NOT_ANNOTATED="${ROOT}/project/metadata/not_annotated_connexin_vertebrate.txt"
NOT_ANNOTATED_WITH_CLADE="${ROOT}/project/metadata/not_annotated_connexin_vertebrate_by_clade.csv"

mkdir -p "${OUT_DIR}" "${LOG_DIR}"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

pids=()
for batch in "${BATCH_DIR}"/*.txt; do
    name=$(basename "${batch}" .txt)
    log "START ${name}"
    (
        python3 "${ROOT}/tools/classify_family_annotation.py" \
            --families connexin \
            --file "${batch}" \
            -o "${OUT_DIR}/status_${name}.csv" \
            --summary "${OUT_DIR}/summary_${name}.csv" \
            >"${LOG_DIR}/${name}.log" 2>&1
    ) &
    pids+=($!)
done

fail=0
for pid in "${pids[@]}"; do
    if ! wait "${pid}"; then
        fail=$((fail + 1))
    fi
done

log "Merging ${#pids[@]} batch outputs (failures=${fail})"
python3 - <<'PY'
import csv
from pathlib import Path

root = Path(__import__("os").environ["ROOT"])
parts = root / "project/metadata/connexin_vertebrate_classify_parts"
catalog = {
    row["organism"]: row["clade"]
    for row in csv.DictReader(
        (root / "project/metadata/species_clade_catalog_connexin_vertebrate.csv").open()
    )
}

status_rows = []
for path in sorted(parts.glob("status_*.csv")):
    with path.open(encoding="utf-8") as handle:
        status_rows.extend(list(csv.DictReader(handle)))

if not status_rows:
    raise SystemExit("No classification results to merge")

status_fields = list(status_rows[0].keys())
merged_status = root / "project/metadata/family_annotation_status_connexin_vertebrate.csv"
with merged_status.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=status_fields)
    writer.writeheader()
    writer.writerows(status_rows)

# wide summary per organism
by_org: dict[str, dict[str, str]] = {}
for row in status_rows:
    org = row["organism"]
    by_org.setdefault(org, {"organism": org})
    by_org[org]["connexin_category"] = row["final_category"]
    by_org[org]["connexin_term"] = row.get("top_matching_term", "")

merged_summary = root / "project/metadata/family_annotation_summary_connexin_vertebrate.csv"
with merged_summary.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(handle, fieldnames=["organism", "connexin_category", "connexin_term"])
    writer.writeheader()
    for org in sorted(by_org):
        writer.writerow(by_org[org])

annotated = "annotated_family_member_found"
not_annotated = [
    org for org, data in sorted(by_org.items())
    if data["connexin_category"] != annotated
]
(root / "project/metadata/not_annotated_connexin_vertebrate.txt").write_text(
    "\n".join(not_annotated) + ("\n" if not_annotated else ""),
    encoding="utf-8",
)

clade_out = root / "project/metadata/not_annotated_connexin_vertebrate_by_clade.csv"
with clade_out.open("w", encoding="utf-8", newline="") as handle:
    writer = csv.DictWriter(
        handle,
        fieldnames=["clade", "organism", "connexin_category", "connexin_term"],
    )
    writer.writeheader()
    for org in not_annotated:
        data = by_org[org]
        writer.writerow(
            {
                "clade": catalog.get(org, ""),
                "organism": org,
                "connexin_category": data["connexin_category"],
                "connexin_term": data["connexin_term"],
            }
        )

print(f"Merged {len(status_rows)} status rows")
print(f"Not annotated connexin: {len(not_annotated)} species")
for org in not_annotated:
    print(f"  - {catalog.get(org,'?')}: {org} ({by_org[org]['connexin_category']})")
PY

log "Done. Lists:"
log "  ${NOT_ANNOTATED}"
log "  ${NOT_ANNOTATED_WITH_CLADE}"
exit "${fail}"
