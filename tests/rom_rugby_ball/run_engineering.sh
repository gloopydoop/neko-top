#!/usr/bin/env bash
set -euo pipefail

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
source "${ROOT}/prepare.env"

CASE_FILE=${CASE_FILE:-"${ROOT}/examples/ROM_rugby_ball/ROM_rugby_ball.case"}
MODE_FILE=${MODE_FILE:-"${ROOT}/results/ROM_rugby_ball/POD_modes/POD_modes0.fld"}
STABLE_MODE_FILE="${ROOT}/results/ROM_rugby_ball.previous_debug_20261117/POD_modes/POD_modes0.fld"
mode_available() {
    local mode_path="$1"
    [ -f "${mode_path}.f00000" ] || \
        [ -f "${mode_path%.fld}.f00000" ]
}
if ! mode_available "${MODE_FILE}" && mode_available "${STABLE_MODE_FILE}"; then
    MODE_FILE="${STABLE_MODE_FILE}"
fi
# The repository launcher clears results/<case> before archiving a run. Keep
# this engineering artifact outside that directory so the launcher cannot
# delete the input basis or the diagnostic outputs during cleanup.
OUT_DIR=${OUT_DIR:-"${ROOT}/results/ROM_rugby_ball_engineering_forward_only_verified"}
HORIZON=${HORIZON:-0.02}
OUTER_STEPS=${OUTER_STEPS:-1}
MPI_RANKS=${MPI_RANKS:-1}
RESET_OUTPUT=${RESET_OUTPUT:-1}

cmake --build "${ROOT}/build" \
    --target ROM_rugby_ball_production_rom_step_diagnostic -j6
ENGINEERING_EXE="${ROOT}/build/tests/rom_rugby_ball/ROM_rugby_ball_production_rom_step_diagnostic"

mkdir -p "${OUT_DIR}"
if [ "${RESET_OUTPUT}" = "1" ]; then
    rm -f "${OUT_DIR}/engineering_summary.csv" \
        "${OUT_DIR}/engineering_rom_status.txt" \
        "${OUT_DIR}/rom_validation.csv" \
        "${OUT_DIR}"/ROM_forward_fields_iter_* \
        "${OUT_DIR}"/ROM_design_iter_* \
        "${OUT_DIR}"/ROM_sensitivity_iter_*
fi
ln -sfn "${ROOT}/box.nmsh" "${OUT_DIR}/box.nmsh"

(
    cd "${OUT_DIR}"
    mpirun -n "${MPI_RANKS}" \
        "${ENGINEERING_EXE}" \
        "${CASE_FILE}" "${MODE_FILE}" file engineering "${HORIZON}" "${OUTER_STEPS}"
)
