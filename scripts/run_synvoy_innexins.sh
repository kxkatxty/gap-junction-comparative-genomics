#!/usr/bin/env bash
# Run SynVoy (Easy Mode) on the deduplicated innexin reference queue.
# Each job searches up to MAX_GENOMES related assemblies for orthologs.
set -uo pipefail
source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
SYNVOY_DIR="${PROJECT_ROOT}/SynVoy"
QUEUE="${PROJECT_ROOT}/project/metadata/innexin_synvoy_queue.tsv"
LOG_DIR="${PROJECT_ROOT}/project/results/synvoy_innexins/logs"
MASTER_LOG="${LOG_DIR}/master.log"
MAX_GENOMES=3  # 3 is faster locally; 5 was reliable but ~2x more download/search work

# Lighter search settings for local runs (saves ~30-40% per job vs defaults).
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
    local outdir="${SYNVOY_DIR}/results/innexin_synvoy/${outdir_name}"
    ls "${outdir}"/synteny_block_*_synteny_plot.html >/dev/null 2>&1
}

# shellcheck disable=SC1091
source "${HOME}/miniconda3/etc/profile.d/conda.sh"
conda activate synvoy_env || { log "FATAL: cannot activate synvoy_env"; exit 1; }
export PATH="${CONDA_PREFIX}/bin:${PATH}"

log "Activated synvoy_env: nextflow $(nextflow -version 2>&1 | grep -i version | head -1)"
cd "${SYNVOY_DIR}" || { log "FATAL: cannot cd to ${SYNVOY_DIR}"; exit 1; }

total=$(($(wc -l < "${QUEUE}") - 1))
idx=0
ok=0
fail=0

tail -n +2 "${QUEUE}" | while IFS=$'\t' read -r accession gene species outdir_name; do
    [ -z "${accession}" ] && continue
    idx=$((idx + 1))
    outdir="results/innexin_synvoy/${outdir_name}"
    run_log="${LOG_DIR}/${outdir_name}.log"

    if job_done "${outdir_name}"; then
        ok=$((ok + 1))
        log "[${idx}/${total}] SKIP  ${gene} (${accession}) already complete -> ${outdir}"
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
        log "[${idx}/${total}] OK    ${gene} (${accession})"
    else
        fail=$((fail + 1))
        log "[${idx}/${total}] FAIL  ${gene} (${accession}) exit=${status} (see ${run_log})"
    fi
done

log "Finished SynVoy innexin queue: ${ok} ok, ${fail} failed. Results under ${SYNVOY_DIR}/results/innexin_synvoy/"
