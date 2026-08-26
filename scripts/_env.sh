# Shared repo paths for scripts/*.sh
# Usage: source "$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)/_env.sh"
_SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
ROOT="$(cd "${_SCRIPT_DIR}/.." && pwd)"
PROJECT_ROOT="$ROOT"
export ROOT PROJECT_ROOT
cd "$ROOT"
export PYTHONPATH="${ROOT}:${ROOT}/tools:${PYTHONPATH:-}"
