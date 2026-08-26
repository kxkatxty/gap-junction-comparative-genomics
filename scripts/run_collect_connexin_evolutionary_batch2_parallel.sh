#!/usr/bin/env bash
# Sample next 5 unannotated-connexin species per evolutionary clade (batch 2).
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
PARTS="${ROOT}/project/metadata/connexin_evolutionary_batch2_parts"
LOG_DIR="${ROOT}/project/results/connexin_evolutionary_collect_batch2/logs"
CATALOG="${ROOT}/project/metadata/species_clade_catalog_connexin_evolutionary_batch2.csv"
QUEUE="${ROOT}/project/metadata/not_annotated_connexin_evolutionary_batch2.txt"
STATUS="${ROOT}/project/metadata/family_annotation_status_connexin_evolutionary_batch2.csv"
EXCLUDE_QUEUE="${ROOT}/project/metadata/not_annotated_connexin_evolutionary.txt"
EXCLUDE_CATALOG="${ROOT}/project/metadata/species_clade_catalog_connexin_evolutionary.csv"

mkdir -p "${PARTS}" "${LOG_DIR}" "${ROOT}/project/results/connexin_evolutionary_collect_batch2"

log() { echo "[$(date '+%H:%M:%S')] $*"; }

CLADES=(
  Chondrichthyes
  Actinopterygii
  Amphibia
  Sauropsida
  Mammalia
)

pids=()
for clade in "${CLADES[@]}"; do
  slug=$(echo "${clade}" | tr '[:upper:]' '[:lower:]')
  log "START ${clade}"
  (
    python3 "${ROOT}/tools/collect_unannotated_connexin_by_clade.py" \
      --clade "${clade}" \
      --per-clade-limit 5 \
      --exclude-txt "${EXCLUDE_QUEUE}" \
      --exclude-catalog "${EXCLUDE_CATALOG}" \
      --output "${PARTS}/catalog_${slug}.csv" \
      --queue "${PARTS}/queue_${slug}.txt" \
      --status "${PARTS}/status_${slug}.csv" \
      >"${LOG_DIR}/${slug}.log" 2>&1
  ) &
  pids+=($!)
done

fail=0
for pid in "${pids[@]}"; do
  if ! wait "${pid}"; then
    fail=$((fail + 1))
  fi
done

log "Merging (${fail} clade job failures)"
python3 - <<'PY'
import csv
from pathlib import Path

root = Path(__import__("os").environ["ROOT"])
parts = root / "project/metadata/connexin_evolutionary_batch2_parts"

catalog_rows = []
queue = []
status_rows = []
for catalog in sorted(parts.glob("catalog_*.csv")):
    text = catalog.read_text(encoding="utf-8").strip()
    if not text:
        continue
    rows = list(csv.DictReader(text.splitlines()))
    catalog_rows.extend(rows)
    queue.extend(row["organism"] for row in rows)

for status in sorted(parts.glob("status_*.csv")):
    text = status.read_text(encoding="utf-8").strip()
    if not text:
        continue
    status_rows.extend(list(csv.DictReader(text.splitlines())))

catalog_out = root / "project/metadata/species_clade_catalog_connexin_evolutionary_batch2.csv"
if catalog_rows:
    fields = list(catalog_rows[0].keys())
    with catalog_out.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields)
        w.writeheader()
        w.writerows(catalog_rows)

queue_out = root / "project/metadata/not_annotated_connexin_evolutionary_batch2.txt"
queue_out.write_text("\n".join(queue) + ("\n" if queue else ""), encoding="utf-8")

status_out = root / "project/metadata/family_annotation_status_connexin_evolutionary_batch2.csv"
if status_rows:
    fields = list(status_rows[0].keys())
    with status_out.open("w", encoding="utf-8", newline="") as handle:
        w = csv.DictWriter(handle, fieldnames=fields)
        w.writeheader()
        w.writerows(status_rows)

print(f"Catalog: {len(catalog_rows)} species across clades")
for clade in sorted({r['clade'] for r in catalog_rows}):
    n = sum(1 for r in catalog_rows if r['clade'] == clade)
    print(f"  {clade}: {n}")
print(f"Queue: {queue_out}")
PY

log "Done -> ${QUEUE}"
exit "${fail}"
