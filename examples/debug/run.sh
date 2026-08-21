#!/bin/bash
set -uo pipefail

if [ -z "${MAIN_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
fi

CASE_FILE=${CASE_FILE:-debug.case}
PYTHON_SCRIPT=${PYTHON_SCRIPT:-"${MAIN_DIR}/scripts/python/pod_state_recover.py"}
NEKO_EXE=${NEKO_BIN:-./neko}
NEKO_RANKS=${NEKO_RANKS:-8}
PY_RANKS=${PY_RANKS:-48}
TOTAL_RANKS=$((NEKO_RANKS + PY_RANKS))
SHORT_TIMEOUT=${DEBUG_TIMEOUT_SHORT:-120}
LONG_TIMEOUT=${DEBUG_TIMEOUT_LONG:-180}
TEST_DIR=${DEBUG_TEST_DIR:-./debug_artifacts}
NO_POD_CASE="${TEST_DIR}/debug_no_pod.case"

mkdir -p "${TEST_DIR}"

source "${MAIN_DIR}/scripts/mpmd_run_helpers.sh"

declare -a STEP_NAMES=()
declare -a STEP_RCS=()

function print_section() {
    printf "\n"
    printf "=%.0s" {1..80}
    printf "\n%s\n" "$1"
    printf "=%.0s" {1..80}
    printf "\n"
}

function record_step_rc() {
    local name=$1
    local rc=$2
    STEP_NAMES+=("${name}")
    STEP_RCS+=("${rc}")
}

function run_step() {
    local name=$1
    shift
    local log_file="${TEST_DIR}/${name}.log"
    local rc

    print_section "STEP ${name}"
    {
        printf 'Command:'
        printf ' %q' "$@"
        printf '\n'
    } > "${log_file}"

    "$@" >> "${log_file}" 2>&1
    rc=$?

    printf 'RC: %s\n' "${rc}" >> "${log_file}"
    cat "${log_file}"
    record_step_rc "${name}" "${rc}"
}

print_section "DEBUG SETUP"

runtime_log="${TEST_DIR}/prepare_python_runtime.log"
{
    printf 'Command: mpmd_prepare_python_runtime %q\n' "${MAIN_DIR}"
    mpmd_prepare_python_runtime "${MAIN_DIR}"
} > "${runtime_log}" 2>&1
runtime_rc=$?
printf 'RC: %s\n' "${runtime_rc}" >> "${runtime_log}"
cat "${runtime_log}"
record_step_rc "prepare_python_runtime" "${runtime_rc}"

env_log="${TEST_DIR}/environment.log"
{
    printf 'Date: %s\n' "$(date --iso-8601=seconds)"
    printf 'PWD: %s\n' "$(pwd)"
    printf 'Host: %s\n' "$(hostname)"
    printf 'Case file: %s\n' "${CASE_FILE}"
    printf 'Neko executable: %s\n' "${NEKO_EXE}"
    printf 'Python script: %s\n' "${PYTHON_SCRIPT}"
    printf 'Neko ranks: %s\n' "${NEKO_RANKS}"
    printf 'Python ranks: %s\n' "${PY_RANKS}"
    printf 'Total ranks: %s\n' "${TOTAL_RANKS}"
    printf '\n--- Module List ---\n'
    module list 2>&1 || true
    printf '\n--- Runtime Environment ---\n'
    mpmd_print_runtime_env
    printf '\n--- SLURM Environment ---\n'
    env | sort | grep '^SLURM_' || true
    printf '\n--- Paths ---\n'
    command -v python3 || true
    command -v python || true
    command -v srun || true
    command -v mpicc || true
    command -v cc || true
    printf '\n--- Directory Listing ---\n'
    ls -la
    printf '\n--- ulimit ---\n'
    ulimit -a
    if [ -n "${SLURM_JOB_ID:-}" ]; then
        printf '\n--- scontrol show job ---\n'
        scontrol show job "${SLURM_JOB_ID}" || true
    fi
} > "${env_log}" 2>&1
cat "${env_log}"
record_step_rc "environment" 0

MPI4PY_SO=$("${PYTHON_BIN}" - <<'PY'
import pathlib
import mpi4py.MPI
print(pathlib.Path(mpi4py.MPI.__file__).resolve())
PY
)

ADIOS2_BINDINGS=$("${PYTHON_BIN}" - <<'PY'
import pathlib
import adios2.bindings
print(pathlib.Path(adios2.bindings.__file__).resolve())
PY
)

cat > "${TEST_DIR}/py_mpi_a.py" <<'EOF'
from mpi4py import MPI
comm = MPI.COMM_WORLD
print(f"A rank {comm.rank} / {comm.size}", flush=True)
EOF

cat > "${TEST_DIR}/py_mpi_b.py" <<'EOF'
from mpi4py import MPI
comm = MPI.COMM_WORLD
print(f"B rank {comm.rank} / {comm.size}", flush=True)
EOF

cat > "${TEST_DIR}/select_gpu.sh" <<'EOF'
#!/bin/bash
export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
exec "$@"
EOF
chmod +x "${TEST_DIR}/select_gpu.sh"

cat > "${TEST_DIR}/select_gpu_delay.sh" <<'EOF'
#!/bin/bash
export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
sleep "${NEKO_STARTUP_DELAY:-20}"
exec "$@"
EOF
chmod +x "${TEST_DIR}/select_gpu_delay.sh"

run_step "prepare_no_pod_case" "${PYTHON_BIN}" - "${CASE_FILE}" "${NO_POD_CASE}" <<'PY'
import json
import pathlib
import sys

src = pathlib.Path(sys.argv[1])
dst = pathlib.Path(sys.argv[2])

with src.open() as f:
    case_data = json.load(f)

case_data["state_recovery"] = {
    "type": "checkpoint",
    "enabled": False,
}

with dst.open("w") as f:
    json.dump(case_data, f, indent=4)
    f.write("\n")

print(dst)
PY

cat > "${TEST_DIR}/mpmd_py_only.conf" <<EOF
0-7 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_a.py
8-55 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

cat > "${TEST_DIR}/mpmd_neko_no_pod_pyhello.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/select_gpu_delay.sh ${NEKO_EXE} ${NO_POD_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

cat > "${TEST_DIR}/mpmd_neko_pyhello.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/select_gpu_delay.sh ${NEKO_EXE} ${CASE_FILE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

cat > "${TEST_DIR}/mpmd_pod_driver.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/select_gpu_delay.sh ${NEKO_EXE} ${CASE_FILE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${PYTHON_SCRIPT} ${CASE_FILE}
EOF

run_step "python_imports" "${PYTHON_BIN}" - <<'PY'
from mpi4py import MPI
import adios2.bindings
print(MPI.Get_library_version())
print(adios2.bindings.__file__)
PY

run_step "mpi4py_ldd" ldd "${MPI4PY_SO}"
run_step "neko_readelf" bash -lc "readelf -d \"${NEKO_EXE}\" | grep -E 'NEEDED|RPATH|RUNPATH'"
run_step "adios2_core_readelf" bash -lc "readelf -d \"${ADIOS2_DIR}/lib64/libadios2_core_mpi.so\" | grep -E 'NEEDED|RPATH|RUNPATH'"
run_step "single_program_python_mpi" timeout --foreground "${SHORT_TIMEOUT}s" \
    srun -n "${TOTAL_RANKS}" --label "${PYTHON_BIN}" "${TEST_DIR}/py_mpi_a.py"
run_step "mpmd_python_only" timeout --foreground "${SHORT_TIMEOUT}s" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog "${TEST_DIR}/mpmd_py_only.conf"
run_step "neko_no_pod_alone" timeout --foreground "${SHORT_TIMEOUT}s" \
    srun -n "${NEKO_RANKS}" --label "${TEST_DIR}/select_gpu.sh" "${NEKO_EXE}" "${NO_POD_CASE}"
run_step "neko_no_pod_pyhello_mpmd" timeout --foreground "${SHORT_TIMEOUT}s" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog "${TEST_DIR}/mpmd_neko_no_pod_pyhello.conf"
run_step "neko_alone" timeout --foreground "${SHORT_TIMEOUT}s" \
    srun -n "${NEKO_RANKS}" --label "${TEST_DIR}/select_gpu.sh" "${NEKO_EXE}" "${CASE_FILE}"
run_step "neko_pyhello_mpmd" timeout --foreground "${SHORT_TIMEOUT}s" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog "${TEST_DIR}/mpmd_neko_pyhello.conf"
run_step "pod_driver_mpmd" timeout --foreground "${LONG_TIMEOUT}s" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog "${TEST_DIR}/mpmd_pod_driver.conf"

summary_file="${TEST_DIR}/summary.txt"
print_section "SUMMARY"
{
    printf '%-28s %s\n' "step" "rc"
    for idx in "${!STEP_NAMES[@]}"; do
        printf '%-28s %s\n' "${STEP_NAMES[$idx]}" "${STEP_RCS[$idx]}"
    done
} > "${summary_file}"
cat "${summary_file}"

printf "\nArtifacts written to %s\n" "${TEST_DIR}"
printf "This debug example always exits successfully so the framework archives the logs.\n"

exit 0
