# Shared env for pannexin/scripts
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${_SCRIPT_DIR}/.." && pwd)"
export ROOT PROJECT_ROOT="$ROOT"
cd "$ROOT"
export PYTHONPATH="${ROOT}:${PYTHONPATH:-}"
