#!/usr/bin/env bash
set -euo pipefail

MINI_MAIN_DIR=${MAIN_DIR:-$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)}
CASE_FILE=${CASE_FILE:-ROM_low_Re_mini.case}
PYTHON_SCRIPT=${PYTHON_SCRIPT:-"${MINI_MAIN_DIR}/scripts/python/pod_state_recover.py"}
NEKO_EXE=${NEKO_BIN:-./neko}
NEKO_RANKS=${NEKO_RANKS:-1}
PY_RANKS=${PY_RANKS:-1}
LOG_FILE=${LOG_FILE:-rom_low_Re_mini_mpmd.log}

if [ -f "${MINI_MAIN_DIR}/prepare.env" ]; then
    # shellcheck disable=SC1091
    source "${MINI_MAIN_DIR}/prepare.env"
fi
# shellcheck disable=SC1091
source "${MINI_MAIN_DIR}/scripts/mpmd_run_helpers.sh"

[ -f "${CASE_FILE}" ] || { echo "Case file not found: ${CASE_FILE}" >&2; exit 1; }
[ -f "${PYTHON_SCRIPT}" ] || { echo "Python driver not found: ${PYTHON_SCRIPT}" >&2; exit 1; }
[ -x "${NEKO_EXE}" ] || { echo "Neko executable not found: ${NEKO_EXE}" >&2; exit 1; }

mpmd_ensure_adios2_python "${MINI_MAIN_DIR}"
mpmd_print_runtime_env
mpmd_launch_shared "${CASE_FILE}" "${PYTHON_SCRIPT}" "${NEKO_EXE}" \
    "${PY_RANKS}" "${NEKO_RANKS}" "${LOG_FILE}"
