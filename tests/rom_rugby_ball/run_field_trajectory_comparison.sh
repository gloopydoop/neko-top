#!/usr/bin/env bash
set -euo pipefail

# Run one real full-order forward/adjoint lifecycle followed by one ROM
# forward/adjoint lifecycle with both paths sampled every 0.5 physical time
# units.  The case is copied and modified in the test output directory so the
# production example remains free of diagnostic-only output settings.

ROOT=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
source "${ROOT}/prepare.env"
# prepare.env intentionally defines installation paths as shell variables for
# the repository launcher.  This test invokes the ROM wrapper as a child
# process, so export the same paths explicitly for its runtime validator.
export ADIOS2_DIR NEKO_DIR JSON_FORTRAN_DIR

OUT_DIR=${OUT_DIR:-"${ROOT}/results/ROM_rugby_ball_field_trajectory"}
MPI_RANKS=${MPI_RANKS:-6}
RESET_OUTPUT=${RESET_OUTPUT:-1}
CASE_OUT="${OUT_DIR}/ROM_rugby_ball_field_trajectory.case"

if [ "${RESET_OUTPUT}" = "1" ]; then
    rm -rf -- "${OUT_DIR}"
fi
mkdir -p -- "${OUT_DIR}"
cp -- "${ROOT}/examples/ROM_rugby_ball/ROM_rugby_ball.case" "${CASE_OUT}"
ln -sfn -- "${ROOT}/box.nmsh" "${OUT_DIR}/box.nmsh"

# Enable physical-time field output for the full-order case and one ROM
# reduced trajectory.  Keep only one optimizer outer step and one inner step:
# this test is a field comparison, not an optimization run.
perl -0pi -e 's/"output_control": "never",/"output_control": "simulationtime",\n            "output_value": 0.5,/s' "${CASE_OUT}"
perl -0pi -e 's/"norm_output_control": "tsteps",/"output_control": "simulationtime",\n            "output_value": 0.5,\n            "norm_output_control": "tsteps",/s' "${CASE_OUT}"
perl -0pi -e 's/"inner_iterations": 30,/"inner_iterations": 1,/; s/("strict_admissibility": false)/$1,\n        "write_trajectory_fields": true,\n        "trajectory_output_interval": 0.5/; s/("type": "rom_mma",\n            "max_iterations":) 30/$1 1/' "${CASE_OUT}"

cmake --build "${ROOT}/build" --target ROM_rugby_ball -j6

export NEKO_BIN="${ROOT}/examples/ROM_rugby_ball/neko"
export PY_SCRIPT="${ROOT}/scripts/python/rom_state_recover.py"
export PY_RANKS="${MPI_RANKS}"
export NEKO_RANKS="${MPI_RANKS}"
export LOG_FILE="mpmd_ROM_rugby_ball_field_trajectory.log"

(
    cd -- "${OUT_DIR}"
    "${ROOT}/examples/ROM_rugby_ball/run.sh" "${CASE_OUT}" "${PY_SCRIPT}"
)

# The volume constraint is not a ROM quantity.  Treat a mismatch between the
# live full-order constraint value/derivative and the values used by MMA as a
# test failure, while leaving the state-dependent objective errors visible in
# rom_validation.csv for engineering inspection.
VALIDATION_CSV="${OUT_DIR}/rom_validation.csv"
if [ ! -s "${VALIDATION_CSV}" ]; then
    echo "Missing ROM validation record: ${VALIDATION_CSV}" >&2
    exit 1
fi
awk -F, 'NR == 2 {
    if (($7 + 0.0) > 1.0e-12 || ($9 + 0.0) > 1.0e-12) {
        printf "Design-only volume constraint reuse failed: value error=%s gradient error=%s\n", $7, $9 > "/dev/stderr"
        exit 1
    }
    printf "Design-only volume constraint check: value error=%s gradient error=%s\n", $7, $9
    exit 0
}' "${VALIDATION_CSV}"

printf 'Field trajectory comparison output: %s\n' "${OUT_DIR}"
printf 'Full-order forward:  %s/field_iter00000_0.*\n' "${OUT_DIR}"
printf 'Full-order adjoint:  %s/adjoint_iter00000_0.*\n' "${OUT_DIR}"
printf 'ROM forward:         %s/ROM_forward_trajectory_iter_*\n' "${OUT_DIR}"
printf 'ROM adjoint:         %s/ROM_adjoint_trajectory_iter_*\n' "${OUT_DIR}"
