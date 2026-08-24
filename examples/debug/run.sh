#!/bin/bash
set -uo pipefail

if [ -z "${MAIN_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
fi

CASE_FILE=${CASE_FILE:-debug.case}
PYTHON_SCRIPT=${PYTHON_SCRIPT:-"${MAIN_DIR}/scripts/python/pod_state_recover.py"}
CURRENT_NEKO_EXE=${NEKO_BIN:-./neko}
NEKO_RANKS=${NEKO_RANKS:-8}
PY_RANKS=${PY_RANKS:-48}
TOTAL_RANKS=$((NEKO_RANKS + PY_RANKS))
DEBUG_USE_GPU=${DEBUG_USE_GPU:-1}
SHORT_TIMEOUT=${DEBUG_TIMEOUT_SHORT:-60}
LONG_TIMEOUT=${DEBUG_TIMEOUT_LONG:-120}
ATTEMPT_TIMEOUT=${DEBUG_TIMEOUT_ATTEMPT:-60}
SMOKE_STEPS=${DEBUG_SMOKE_STEPS:-1}
DOCS_CPUS_PER_TASK=${DEBUG_CPUS_PER_TASK:-7}
DOCS_GPU_MAP=${DEBUG_LUMI_GPU_MAP:-4,5,2,3,6,7,0,1}
DOCS_CPU_BIND=${DEBUG_LUMI_CPU_BIND:-cores}
DOCS_CPU_MASK=${DEBUG_LUMI_CPU_MASK:-mask_cpu:7e000000000000,7e00000000000000,7e0000,7e000000,7e,7e00,7e00000000,7e0000000000}
SMALL_PY_RANKS=${DEBUG_SMALL_PY_RANKS:-8}
SMALL_TOTAL_RANKS=$((NEKO_RANKS + SMALL_PY_RANKS))
TEST_DIR=${DEBUG_TEST_DIR:-./debug_artifacts}
SMOKE_CASE="${TEST_DIR}/debug_smoke.case"
NO_POD_CASE="${TEST_DIR}/debug_smoke_no_pod.case"

mkdir -p "${TEST_DIR}"

source "${MAIN_DIR}/scripts/mpmd_run_helpers.sh"
RUNTIME_ENV_FILE=${MPMD_RUNTIME_ENV_FILE:-"${MAIN_DIR}/build/mpmd_runtime.env"}

if [ "${DEBUG_USE_GPU}" = "1" ]; then
    STANDALONE_DOCS_ARGS=(
        --exact -n "${NEKO_RANKS}" --cpus-per-task="${DOCS_CPUS_PER_TASK}"
        --cpu-bind="${DOCS_CPU_BIND}" --gpu-bind="map:${DOCS_GPU_MAP}"
        --gres-flags=allow-task-sharing --label
    )
    STANDALONE_WRAPPER_ARGS=(
        --exact -n "${NEKO_RANKS}" --cpus-per-task="${DOCS_CPUS_PER_TASK}"
        --cpu-bind="${DOCS_CPU_BIND}" --gres-flags=allow-task-sharing --label
    )
else
    STANDALONE_DOCS_ARGS=(
        --exact -n "${NEKO_RANKS}" --cpus-per-task="${DOCS_CPUS_PER_TASK}"
        --cpu-bind="${DOCS_CPU_BIND}" --label
    )
    STANDALONE_WRAPPER_ARGS=("${STANDALONE_DOCS_ARGS[@]}")
fi

declare -a STEP_NAMES=()
declare -a STEP_RCS=()
declare -A STEP_RC_MAP=()

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
    STEP_RC_MAP["${name}"]="${rc}"
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
    return 0
}

function reset_neko_outputs() {
    rm -rf ./checkpoints ./design ./forward_fields ./adjoint_fields ./sensitivity
    rm -f ./adjoint_norm.csv ./optimization_data.csv ./joblimit*.chkp
    rm -f ./core ./core.*
}

function run_neko_step() {
    local name=$1
    local timeout_seconds=$2
    shift 2

    reset_neko_outputs
    run_step "${name}" timeout --foreground "${timeout_seconds}s" "$@"
}

function skip_step() {
    local name=$1
    local reason=$2
    local log_file="${TEST_DIR}/${name}.log"

    print_section "STEP ${name}"
    printf 'SKIP: %s\n' "${reason}" > "${log_file}"
    cat "${log_file}"
    record_step_rc "${name}" "SKIP"
}

function finish_debug() {
    local summary_file="${TEST_DIR}/summary.txt"

    print_section "SUMMARY"
    {
        printf '%-32s %s\n' "step" "rc"
        for idx in "${!STEP_NAMES[@]}"; do
            printf '%-32s %s\n' "${STEP_NAMES[$idx]}" "${STEP_RCS[$idx]}"
        done
    } > "${summary_file}"
    cat "${summary_file}"

    printf "\nArtifacts written to %s\n" "${TEST_DIR}"
    printf "This reduced debug example exits successfully so the framework archives the logs.\n"
    exit 0
}

function resolve_optional_binary() {
    local candidate

    for candidate in "$@"; do
        [ -n "${candidate}" ] || continue
        if [ -x "${candidate}" ]; then
            realpath "${candidate}"
            return 0
        fi
    done

    return 1
}

if [ -x "${CURRENT_NEKO_EXE}" ]; then
    CURRENT_NEKO_EXE=$(realpath "${CURRENT_NEKO_EXE}")
fi
WORKING_NEKO_EXE=$(resolve_optional_binary \
    "${DEBUG_WORKING_NEKO:-}" \
    "./working_neko" \
    "/scratch/nobis/debug/neko" 2>/dev/null || true)
NODEVICE_NEKO_EXE=$(resolve_optional_binary \
    "${DEBUG_NODEVICE_NEKO:-}" \
    "./nodevice_neko" 2>/dev/null || true)
BAD_NEKO_EXE=$(resolve_optional_binary \
    "${DEBUG_BAD_NEKO:-}" \
    "./bad_neko" \
    "/scratch/nobis/debug/bad_neko" 2>/dev/null || true)

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
    printf 'Current Neko: %s\n' "${CURRENT_NEKO_EXE}"
    printf 'Working Neko: %s\n' "${WORKING_NEKO_EXE:-<unavailable>}"
    printf 'No-device Neko: %s\n' "${NODEVICE_NEKO_EXE:-<unavailable>}"
    printf 'Bad Neko: %s\n' "${BAD_NEKO_EXE:-<unavailable>}"
    printf 'Python script: %s\n' "${PYTHON_SCRIPT}"
    printf 'Neko ranks: %s\n' "${NEKO_RANKS}"
    printf 'Python ranks: %s\n' "${PY_RANKS}"
    printf 'Total ranks: %s\n' "${TOTAL_RANKS}"
    printf 'Debug use GPU: %s\n' "${DEBUG_USE_GPU}"
    printf 'Smoke steps: %s\n' "${SMOKE_STEPS}"
    printf 'LUMI GPU map: %s\n' "${DOCS_GPU_MAP}"
    printf 'LUMI cpu-bind: %s\n' "${DOCS_CPU_BIND}"
    printf 'LUMI cpus-per-task: %s\n' "${DOCS_CPUS_PER_TASK}"
    printf 'Runtime env: %s\n' "${RUNTIME_ENV_FILE}"
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

cat > "${TEST_DIR}/py_mpi_b.py" <<'EOF'
from mpi4py import MPI
comm = MPI.COMM_WORLD
print(f"B rank {comm.rank} / {comm.size}", flush=True)
EOF

cat > "${TEST_DIR}/py_split_peer.py" <<'EOF'
#!/usr/bin/env python3
import os
import time
from mpi4py import MPI

world = MPI.COMM_WORLD
color = int(os.environ["NEKO_COMM_ID"])
local = world.Split(color, world.rank)

print(
    "split_peer "
    f"world_rank={world.rank} world_size={world.size} "
    f"local_rank={local.rank} local_size={local.size} color={color}",
    flush=True,
)

time.sleep(int(os.environ.get("SPLIT_PEER_SLEEP_SECONDS", "20")))
EOF
chmod +x "${TEST_DIR}/py_split_peer.py"

cat > "${TEST_DIR}/helper_sleep.sh" <<'EOF'
#!/bin/bash
set -eu
echo "helper rank ${SLURM_PROCID:-?} / local ${SLURM_LOCALID:-?}" >&2
sleep "${HELPER_SLEEP_SECONDS:-30}"
EOF
chmod +x "${TEST_DIR}/helper_sleep.sh"

cat > "${TEST_DIR}/run_neko_env.sh" <<'EOF'
#!/bin/bash
export MPICH_DMAPP_APP_IS_WORLD=${MPICH_DMAPP_APP_IS_WORLD:-1}
export MPICH_GPU_SUPPORT_ENABLED=${MPICH_GPU_SUPPORT_ENABLED:-1}
export NEKO_GS_STRTGY=${NEKO_GS_STRTGY:-3}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-2}
export OMP_WAIT_POLICY=${OMP_WAIT_POLICY:-PASSIVE}
exec "$@"
EOF
chmod +x "${TEST_DIR}/run_neko_env.sh"

cat > "${TEST_DIR}/run_neko_mpmd_env.sh" <<'EOF'
#!/bin/bash
export MPICH_DMAPP_APP_IS_WORLD=${MPICH_DMAPP_APP_IS_WORLD:-1}
export MPICH_GPU_SUPPORT_ENABLED=${MPICH_GPU_SUPPORT_ENABLED:-1}
export NEKO_GS_STRTGY=${NEKO_GS_STRTGY:-3}
export NEKO_GS_COMM=${NEKO_GS_COMM:-MPI}
export NEKO_DISABLE_DEVICE_MPI=${NEKO_DISABLE_DEVICE_MPI:-1}
export OMP_NUM_THREADS=${OMP_NUM_THREADS:-1}
export OMP_WAIT_POLICY=${OMP_WAIT_POLICY:-PASSIVE}
exec "$@"
EOF
chmod +x "${TEST_DIR}/run_neko_mpmd_env.sh"

cat > "${TEST_DIR}/run_python_host_env.sh" <<'EOF'
#!/bin/bash
unset ROCR_VISIBLE_DEVICES
export MPICH_GPU_SUPPORT_ENABLED=0
exec "$@"
EOF
chmod +x "${TEST_DIR}/run_python_host_env.sh"

cat > "${TEST_DIR}/select_gpu_preserve.sh" <<'EOF'
#!/bin/bash
if [ "${DEBUG_USE_GPU:-1}" = "1" ] && [ -z "${ROCR_VISIBLE_DEVICES:-}" ]; then
    export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
fi
exec "$@"
EOF
chmod +x "${TEST_DIR}/select_gpu_preserve.sh"

run_step "prepare_smoke_cases" "${PYTHON_BIN}" - "${CASE_FILE}" "${SMOKE_CASE}" "${NO_POD_CASE}" "${SMOKE_STEPS}" <<'PY'
import json
import pathlib
import sys

src = pathlib.Path(sys.argv[1])
pod_dst = pathlib.Path(sys.argv[2])
no_pod_dst = pathlib.Path(sys.argv[3])
smoke_steps = max(int(sys.argv[4]), 1)

with src.open() as f:
    base = json.load(f)

dt = float(base["case"]["time"]["timestep"])
end_time = dt * smoke_steps

def make_common(data):
    data["case"]["time"]["end_time"] = end_time
    solver = data.get("optimization", {}).get("solver", {})
    solver["max_iterations"] = 1
    for objective in data.get("optimization", {}).get("objectives", []):
        if "start_time" in objective:
            objective["start_time"] = 0.0
    return data

pod_case = make_common(json.loads(json.dumps(base)))
state_recovery = pod_case.setdefault("state_recovery", {})
state_recovery["enabled"] = True
state_recovery["type"] = "pod"
state_recovery["i_stream"] = 1
state_recovery["batch_size"] = 1
state_recovery["n_memory"] = max(int(state_recovery.get("n_memory", 2)), 2)
state_recovery["n_modes"] = min(int(state_recovery.get("n_modes", 1)), 4)

no_pod_case = make_common(json.loads(json.dumps(base)))
no_pod_case["state_recovery"] = {
    "type": "checkpoint",
    "enabled": False,
}

for path, payload in ((pod_dst, pod_case), (no_pod_dst, no_pod_case)):
    with path.open("w") as f:
        json.dump(payload, f, indent=4)
        f.write("\n")

print(f"pod_case={pod_dst}")
print(f"no_pod_case={no_pod_dst}")
print(f"smoke_end_time={end_time}")
PY

run_step "verify_current_neko" test -x "${CURRENT_NEKO_EXE}"
if [ "${STEP_RC_MAP[verify_current_neko]:-1}" != "0" ]; then
    finish_debug
fi

cat > "${TEST_DIR}/mpmd_current_no_pod.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_mpmd_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${CURRENT_NEKO_EXE} ${NO_POD_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

cat > "${TEST_DIR}/mpmd_current_no_pod_helpers.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_mpmd_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${CURRENT_NEKO_EXE} ${NO_POD_CASE}
8-55 ${TEST_DIR}/helper_sleep.sh
EOF

cat > "${TEST_DIR}/mpmd_current_no_pod_splitpeer.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_mpmd_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${CURRENT_NEKO_EXE} ${NO_POD_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${TEST_DIR}/py_split_peer.py
EOF

cat > "${TEST_DIR}/mpmd_current_no_pod_hostpy.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_mpmd_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${CURRENT_NEKO_EXE} ${NO_POD_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${TEST_DIR}/run_python_host_env.sh ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

run_neko_step "current_no_pod_alone_docs" "${SHORT_TIMEOUT}" \
    srun "${STANDALONE_DOCS_ARGS[@]}" \
    "${TEST_DIR}/run_neko_env.sh" "${CURRENT_NEKO_EXE}" "${NO_POD_CASE}"
if [ "${STEP_RC_MAP[current_no_pod_alone_docs]:-1}" != "0" ]; then
    skip_step "current_no_pod_mpmd_helpers" \
        "Skipped because standalone current Neko failed."
    skip_step "current_no_pod_mpmd_splitpeer" \
        "Skipped because standalone current Neko failed."
    skip_step "current_no_pod_mpmd_hostpy" \
        "Skipped because standalone current Neko failed."
    skip_step "current_no_pod_mpmd" \
        "Skipped because standalone current Neko failed."
    finish_debug
fi

run_neko_step "current_no_pod_mpmd_helpers" "${ATTEMPT_TIMEOUT}" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod_helpers.conf"
if [ "${STEP_RC_MAP[current_no_pod_mpmd_helpers]:-1}" = "143" ]; then
    skip_step "current_no_pod_mpmd_splitpeer" \
        "Skipped because helper MPMD hit the PMI assert."
    skip_step "current_no_pod_mpmd_hostpy" \
        "Skipped because helper MPMD hit the PMI assert."
    skip_step "current_no_pod_mpmd" \
        "Skipped because helper MPMD hit the PMI assert."
    finish_debug
fi

run_neko_step "current_no_pod_mpmd_splitpeer" "${ATTEMPT_TIMEOUT}" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod_splitpeer.conf"
if [ "${STEP_RC_MAP[current_no_pod_mpmd_splitpeer]:-1}" != "0" ]; then
    skip_step "current_no_pod_mpmd_hostpy" \
        "Skipped because split-peer MPMD already failed."
    skip_step "current_no_pod_mpmd" \
        "Skipped because split-peer MPMD already failed."
    finish_debug
fi

run_neko_step "current_no_pod_mpmd_hostpy" "${ATTEMPT_TIMEOUT}" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod_hostpy.conf"
if [ "${STEP_RC_MAP[current_no_pod_mpmd_hostpy]:-1}" != "0" ]; then
    skip_step "current_no_pod_mpmd" \
        "Skipped because host-Python MPMD already failed."
    finish_debug
fi

run_neko_step "current_no_pod_mpmd" "${LONG_TIMEOUT}" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod.conf"
finish_debug
