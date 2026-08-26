#!/usr/bin/env bash
# Run SynVoy Easy Mode on a TSV queue (accession, gene, species, outdir_name).
# Usage:
#   ./scripts/run_synvoy_queue.sh project/metadata/connexin_synvoy_queue.tsv connexin_synvoy
#   ./scripts/run_synvoy_queue.sh project/metadata/innexin_synvoy_queue_insects.tsv innexin_synvoy
#   ./scripts/run_synvoy_queue.sh project/metadata/innexin_synvoy_queue_nematodes.tsv innexin_synvoy
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
SYNVOY_DIR="${PROJECT_ROOT}/SynVoy"
QUEUE_ARG="${1:?Usage: $0 <queue.tsv> [results_subdir]}"
# Resolve queue path before cd'ing into SynVoy
if [[ "${QUEUE_ARG}" = /* ]]; then
    QUEUE="${QUEUE_ARG}"
else
    QUEUE="${PROJECT_ROOT}/${QUEUE_ARG}"
fi
if [[ ! -f "${QUEUE}" ]]; then
    echo "FATAL: queue file not found: ${QUEUE}" >&2
    exit 1
fi
RESULTS_SUBDIR="${2:-innexin_synvoy}"
LOG_DIR="${PROJECT_ROOT}/project/results/synvoy_runs/logs/${RESULTS_SUBDIR}"
MASTER_LOG="${LOG_DIR}/master.log"
MAX_GENOMES="${MAX_GENOMES:-3}"

SYNVOY_FAST_ARGS=(
    --n_flanking_genes 5
    --adaptive_max_regions 3
    --enable_smith_waterman false
    --exon_level_search false
)

mkdir -p "${LOG_DIR}"
log() { echo "[$(date '+%Y-%m-%d %H:%M:%S')] $*" | tee -a "${MASTER_LOG}"; }

job_done() {
    local outdir_name="$1"
    local outdir="${SYNVOY_DIR}/results/${RESULTS_SUBDIR}/${outdir_name}"
    ls "${outdir}"/synteny_block_*_synteny_plot.html >/dev/null 2>&1
}

# shellcheck disable=SC1091
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate synvoy_env || { log "FATAL: cannot activate synvoy_env"; exit 1; }
export PATH="${CONDA_PREFIX}/bin:${PATH}"
export NXF_DISABLE_CHECK_LATEST="${NXF_DISABLE_CHECK_LATEST:-true}"
export NXF_DISABLE_CHECK_LATEST="${NXF_DISABLE_CHECK_LATEST:-true}"

cd "${SYNVOY_DIR}" || { log "FATAL: cannot cd to SynVoy"; exit 1; }
log "Queue=${QUEUE} results=results/${RESULTS_SUBDIR} max_genomes=${MAX_GENOMES}"

total=$(($(wc -l < "${QUEUE}") - 1))
idx=0
ok=0
fail=0

while IFS=$'\t' read -r accession gene species outdir_name; do
    [ "${accession}" = "accession" ] && continue
    [ -z "${accession}" ] && continue
    idx=$((idx + 1))
    outdir="results/${RESULTS_SUBDIR}/${outdir_name}"
    run_log="${LOG_DIR}/${outdir_name}.log"

    if job_done "${outdir_name}"; then
        ok=$((ok + 1))
        log "[${idx}/${total}] SKIP  ${gene} (${accession}) already complete"
        continue
    fi

    log "[${idx}/${total}] START ${gene} (${accession}, ${species}) -> ${outdir}"
    nextflow run main.nf \
        -profile standard \
        --mode easy \
        --query_id "${accession}" \
        --max_genomes "${MAX_GENOMES}" \
        --outdir "${outdir}" \
        --auto_params false \
        --multi_profile false \
        "${SYNVOY_FAST_ARGS[@]}" \
        -resume \
        >"${run_log}" 2>&1
    status=$?
    if [ "${status}" -eq 0 ]; then
        ok=$((ok + 1))
        log "[${idx}/${total}] OK    ${gene}"
    else
        fail=$((fail + 1))
        log "[${idx}/${total}] FAIL  ${gene} exit=${status} (see ${run_log})"
    fi
done < "${QUEUE}"

log "Finished: ${ok} ok, ${fail} failed under results/${RESULTS_SUBDIR}/"
