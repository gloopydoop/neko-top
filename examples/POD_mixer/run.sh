#!/bin/bash
set -euo pipefail

if [ -z "${MAIN_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
fi

source "${MAIN_DIR}/scripts/mpmd_run_helpers.sh"

CASE_DIR=$(pwd)
CASE_FILE=$(mpmd_resolve_case_file "${CASE_DIR}" "${1:-}")
PY_SCRIPT=${2:-"${MAIN_DIR}/scripts/python/pod_state_recover.py"}

LOCAL_NEKO_RANKS=${NEKO_RANKS:-6}
LOCAL_PY_RANKS=${PY_RANKS:-6}
NEKO_RANKS_PER_NODE=${NEKO_RANKS_PER_NODE:-8}
PY_RANKS_PER_NODE=${PY_RANKS_PER_NODE:-48}

if [ -x "${CASE_DIR}/prepare.sh" ] && [ ! -f "${CASE_DIR}/box.nmsh" ]; then
    ./prepare.sh
fi

if [ -n "${NEKO_BIN:-}" ]; then
    NEKO_EXE="${NEKO_BIN}"
elif [ -x "./neko" ]; then
    NEKO_EXE="./neko"
elif [ -x "../neko" ]; then
    NEKO_EXE="../neko"
else
    echo "Error: no Neko executable found." >&2
    exit 1
fi

CASE_NAME="$(basename "${CASE_FILE}")"
CASE_BASE="${CASE_NAME%.case}"
LOG_FILE=${LOG_FILE:-"mpmd_${CASE_BASE}.log"}

mpmd_prepare_python_runtime "${MAIN_DIR}"
mpmd_print_runtime_env

if [ -n "${SLURM_JOB_ID:-}" ] && [ -n "${POD_MIXER_MPMD_CONF_GENERATOR:-}" ]; then
    if [ ! -x "${POD_MIXER_MPMD_CONF_GENERATOR}" ]; then
        echo "Error: POD_mixer MPMD generator is not executable:" >&2
        echo "  ${POD_MIXER_MPMD_CONF_GENERATOR}" >&2
        exit 1
    fi

    "${POD_MIXER_MPMD_CONF_GENERATOR}" \
        --nodes "${SLURM_NNODES:-1}" \
        --case-file "${CASE_FILE}" \
        --neko-exe "${NEKO_EXE}" \
        --python-bin "${PYTHON_BIN}" \
        --python-script "${PY_SCRIPT}" \
        --neko-ranks-per-node "${NEKO_RANKS_PER_NODE}" \
        --python-ranks-per-node "${PY_RANKS_PER_NODE}" \
        --output "./mpmd.conf" \
        --select-gpu "./select_gpu"

    trap 'rm -f ./select_gpu ./mpmd.conf' EXIT

    read -r -a srun_extra <<< "${SRUN_MPMD_FLAGS:---distribution=block:block}"

    total_neko_ranks=$((NEKO_RANKS_PER_NODE * ${SLURM_NNODES:-1}))
    total_py_ranks=$((PY_RANKS_PER_NODE * ${SLURM_NNODES:-1}))

    echo "Launching ${total_neko_ranks} Neko GPU ranks and ${total_py_ranks} Python ranks"
    echo "Case:   ${CASE_FILE}"
    echo "Neko:   ${NEKO_EXE}"
    echo "Python: ${PYTHON_BIN} ${PY_SCRIPT}"
    echo "Config: ./mpmd.conf"
    echo "Output: ${LOG_FILE}"

    if ! srun --unbuffered "${srun_extra[@]}" --multi-prog ./mpmd.conf \
        > "${LOG_FILE}" 2>&1; then
        echo "Error: shared MPMD launch failed. See ${LOG_FILE}." >&2
        exit 1
    fi

    exit 0
fi

mpmd_launch_shared "${CASE_FILE}" "${PY_SCRIPT}" "${NEKO_EXE}" \
    "${LOCAL_PY_RANKS}" "${LOCAL_NEKO_RANKS}" "${LOG_FILE}"
