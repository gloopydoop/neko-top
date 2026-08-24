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

cat > "${TEST_DIR}/helper_sleep.sh" <<'EOF'
#!/bin/bash
set -eu
echo "helper rank ${SLURM_PROCID:-?} / local ${SLURM_LOCALID:-?}" >&2
sleep "${HELPER_SLEEP_SECONDS:-30}"
EOF
chmod +x "${TEST_DIR}/helper_sleep.sh"

cat > "${TEST_DIR}/binary_report.sh" <<'EOF'
#!/bin/bash
set -eu

binary=$1

file "${binary}"
ls -lh "${binary}"
sha256sum "${binary}"
strings "${binary}" | grep -m1 '(build:' || true
readelf -d "${binary}" | grep -E 'NEEDED|RPATH|RUNPATH'
EOF
chmod +x "${TEST_DIR}/binary_report.sh"

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

cat > "${TEST_DIR}/select_gpu_force.sh" <<'EOF'
#!/bin/bash
export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
exec "$@"
EOF
chmod +x "${TEST_DIR}/select_gpu_force.sh"

cat > "${TEST_DIR}/select_gpu_preserve.sh" <<'EOF'
#!/bin/bash
if [ -z "${ROCR_VISIBLE_DEVICES:-}" ]; then
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

cat > "${TEST_DIR}/mpmd_py_only.conf" <<EOF
0-7 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_a.py
8-55 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

cat > "${TEST_DIR}/mpmd_current_no_pod.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_mpmd_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${CURRENT_NEKO_EXE} ${NO_POD_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

cat > "${TEST_DIR}/mpmd_current_pod.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_mpmd_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${CURRENT_NEKO_EXE} ${SMOKE_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${PYTHON_SCRIPT} ${SMOKE_CASE}
EOF

cat > "${TEST_DIR}/mpmd_current_no_pod_helpers.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_mpmd_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${CURRENT_NEKO_EXE} ${NO_POD_CASE}
8-55 ${TEST_DIR}/helper_sleep.sh
EOF

cat > "${TEST_DIR}/mpmd_current_no_pod_small.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_mpmd_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${CURRENT_NEKO_EXE} ${NO_POD_CASE}
8-$((SMALL_TOTAL_RANKS - 1)) /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

cat > "${TEST_DIR}/mpmd_current_no_pod_hostpy.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_mpmd_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${CURRENT_NEKO_EXE} ${NO_POD_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${TEST_DIR}/run_python_host_env.sh ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

if [ -n "${WORKING_NEKO_EXE:-}" ]; then
    cat > "${TEST_DIR}/mpmd_working_no_pod.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${WORKING_NEKO_EXE} ${NO_POD_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

    cat > "${TEST_DIR}/mpmd_working_pod.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${WORKING_NEKO_EXE} ${SMOKE_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${PYTHON_SCRIPT} ${SMOKE_CASE}
EOF
fi

if [ -n "${NODEVICE_NEKO_EXE:-}" ]; then
    cat > "${TEST_DIR}/mpmd_nodevice_no_pod.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${NODEVICE_NEKO_EXE} ${NO_POD_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${TEST_DIR}/py_mpi_b.py
EOF

    cat > "${TEST_DIR}/mpmd_nodevice_pod.conf" <<EOF
0-7 /usr/bin/env NEKO_COMM_ID=0 NEKO_CTRL_PEER_ROOT=${NEKO_RANKS} ${TEST_DIR}/run_neko_env.sh ${TEST_DIR}/select_gpu_preserve.sh ${NODEVICE_NEKO_EXE} ${SMOKE_CASE}
8-55 /usr/bin/env NEKO_COMM_ID=1 NEKO_CTRL_PEER_ROOT=0 ${PYTHON_BIN} ${PYTHON_SCRIPT} ${SMOKE_CASE}
EOF
fi

run_step "python_imports" "${PYTHON_BIN}" - <<'PY'
from mpi4py import MPI
import adios2.bindings
print(MPI.Get_library_version())
print(adios2.bindings.__file__)
PY

run_step "srun_mpi_list" srun --mpi=list
run_step "mpi4py_ldd" ldd "${MPI4PY_SO}"
run_step "binary_current" "${TEST_DIR}/binary_report.sh" "${CURRENT_NEKO_EXE}"

if [ -n "${WORKING_NEKO_EXE:-}" ]; then
    run_step "binary_working" "${TEST_DIR}/binary_report.sh" "${WORKING_NEKO_EXE}"
    run_step "binary_compare_current_working" "${PYTHON_BIN}" - "${CURRENT_NEKO_EXE}" "${WORKING_NEKO_EXE}" <<'PY'
import pathlib
import re
import subprocess
import sys

pattern = re.compile(r"Shared library: \[(.*?)\]")
rpath_pattern = re.compile(r"Library (?:rpath|runpath): \[(.*?)\]")

def inspect(path):
    out = subprocess.check_output(["readelf", "-d", path], text=True)
    return {
        "needed": sorted(set(pattern.findall(out))),
        "rpath": rpath_pattern.findall(out),
    }

left = inspect(sys.argv[1])
right = inspect(sys.argv[2])

print(f"current_only={sorted(set(left['needed']) - set(right['needed']))}")
print(f"working_only={sorted(set(right['needed']) - set(left['needed']))}")
print(f"current_rpath={left['rpath']}")
print(f"working_rpath={right['rpath']}")
PY
else
    skip_step "binary_working" "No working control binary available."
    skip_step "binary_compare_current_working" \
        "No working control binary available."
fi

if [ -n "${BAD_NEKO_EXE:-}" ]; then
    run_step "binary_bad" "${TEST_DIR}/binary_report.sh" "${BAD_NEKO_EXE}"
else
    skip_step "binary_bad" "No bad control binary available."
fi

if [ -n "${NODEVICE_NEKO_EXE:-}" ]; then
    run_step "binary_nodevice" "${TEST_DIR}/binary_report.sh" "${NODEVICE_NEKO_EXE}"
else
    skip_step "binary_nodevice" "No no-device-MPI binary available."
fi

run_step "mpmd_python_only" timeout --foreground "${SHORT_TIMEOUT}s" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_py_only.conf"

run_neko_step "current_no_pod_alone_wrapper" "${SHORT_TIMEOUT}" \
    srun -n "${NEKO_RANKS}" --label \
    "${TEST_DIR}/run_neko_env.sh" "${TEST_DIR}/select_gpu_force.sh" \
    "${CURRENT_NEKO_EXE}" "${NO_POD_CASE}"

run_neko_step "current_no_pod_alone_docs" "${SHORT_TIMEOUT}" \
    srun --exact -n "${NEKO_RANKS}" --cpus-per-task="${DOCS_CPUS_PER_TASK}" \
    --cpu-bind="${DOCS_CPU_BIND}" --gpu-bind="map:${DOCS_GPU_MAP}" \
    --gres-flags=allow-task-sharing --label \
    "${TEST_DIR}/run_neko_env.sh" "${CURRENT_NEKO_EXE}" "${NO_POD_CASE}"

run_neko_step "current_no_pod_alone_docs_preserve" "${SHORT_TIMEOUT}" \
    srun --exact -n "${NEKO_RANKS}" --cpus-per-task="${DOCS_CPUS_PER_TASK}" \
    --cpu-bind="${DOCS_CPU_BIND}" --gpu-bind="map:${DOCS_GPU_MAP}" \
    --gres-flags=allow-task-sharing --label \
    "${TEST_DIR}/run_neko_env.sh" "${TEST_DIR}/select_gpu_preserve.sh" \
    "${CURRENT_NEKO_EXE}" "${NO_POD_CASE}"

run_neko_step "current_no_pod_alone_docs_wrapper" "${SHORT_TIMEOUT}" \
    srun --exact -n "${NEKO_RANKS}" --cpus-per-task="${DOCS_CPUS_PER_TASK}" \
    --cpu-bind="${DOCS_CPU_BIND}" --gres-flags=allow-task-sharing --label \
    "${TEST_DIR}/run_neko_env.sh" "${TEST_DIR}/select_gpu_force.sh" \
    "${CURRENT_NEKO_EXE}" "${NO_POD_CASE}"

run_neko_step "current_no_pod_alone_mask_wrapper" "${SHORT_TIMEOUT}" \
    srun --exact -n "${NEKO_RANKS}" --cpus-per-task="${DOCS_CPUS_PER_TASK}" \
    --cpu-bind="${DOCS_CPU_MASK}" --gres-flags=allow-task-sharing --label \
    "${TEST_DIR}/run_neko_env.sh" "${TEST_DIR}/select_gpu_force.sh" \
    "${CURRENT_NEKO_EXE}" "${NO_POD_CASE}"

if [ -n "${WORKING_NEKO_EXE:-}" ]; then
    run_neko_step "working_no_pod_alone_docs" "${SHORT_TIMEOUT}" \
        srun --exact -n "${NEKO_RANKS}" --cpus-per-task="${DOCS_CPUS_PER_TASK}" \
        --cpu-bind="${DOCS_CPU_BIND}" --gpu-bind="map:${DOCS_GPU_MAP}" \
        --gres-flags=allow-task-sharing --label \
        "${TEST_DIR}/run_neko_env.sh" "${WORKING_NEKO_EXE}" "${NO_POD_CASE}"
else
    skip_step "working_no_pod_alone_docs" "No working control binary available."
fi

if [ -n "${NODEVICE_NEKO_EXE:-}" ]; then
    run_neko_step "nodevice_no_pod_alone_docs" "${SHORT_TIMEOUT}" \
        srun --exact -n "${NEKO_RANKS}" --cpus-per-task="${DOCS_CPUS_PER_TASK}" \
        --cpu-bind="${DOCS_CPU_BIND}" --gpu-bind="map:${DOCS_GPU_MAP}" \
        --gres-flags=allow-task-sharing --label \
        "${TEST_DIR}/run_neko_env.sh" "${NODEVICE_NEKO_EXE}" "${NO_POD_CASE}"
else
    skip_step "nodevice_no_pod_alone_docs" "No no-device-MPI binary available."
fi

run_neko_step "current_no_pod_mpmd_helpers" "${ATTEMPT_TIMEOUT}" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod_helpers.conf"

run_neko_step "current_no_pod_mpmd_helpers_shasta" "${ATTEMPT_TIMEOUT}" \
    srun --mpi=cray_shasta -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod_helpers.conf"

run_neko_step "current_no_pod_mpmd_helpers_pmi2" "${ATTEMPT_TIMEOUT}" \
    srun --mpi=pmi2 -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod_helpers.conf"

if [ "${STEP_RC_MAP[current_no_pod_mpmd_helpers]:-1}" != "0" ]; then
    run_neko_step "current_no_pod_mpmd_helpers_env_display" "${ATTEMPT_TIMEOUT}" \
        env MPICH_ENV_DISPLAY=1 MPICH_VERSION_DISPLAY=1 MPICH_CPUMASK_DISPLAY=1 \
        MPICH_ABORT_ON_ERROR=1 \
        srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
        "${TEST_DIR}/mpmd_current_no_pod_helpers.conf"
else
    skip_step "current_no_pod_mpmd_helpers_env_display" \
        "Skipped because current_no_pod_mpmd_helpers already passed."
fi

run_neko_step "current_no_pod_mpmd_small" "${ATTEMPT_TIMEOUT}" \
    srun --exact -n "${SMALL_TOTAL_RANKS}" --label \
    --unbuffered --multi-prog "${TEST_DIR}/mpmd_current_no_pod_small.conf"

run_neko_step "current_no_pod_mpmd_hostpy" "${ATTEMPT_TIMEOUT}" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod_hostpy.conf"

run_neko_step "current_no_pod_mpmd_hostpy_shasta" "${ATTEMPT_TIMEOUT}" \
    srun --mpi=cray_shasta -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod_hostpy.conf"

run_neko_step "current_no_pod_mpmd_hostpy_pmi2" "${ATTEMPT_TIMEOUT}" \
    srun --mpi=pmi2 -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod_hostpy.conf"

run_neko_step "current_no_pod_mpmd" "${LONG_TIMEOUT}" \
    srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
    "${TEST_DIR}/mpmd_current_no_pod.conf"

if [ -n "${WORKING_NEKO_EXE:-}" ]; then
    run_neko_step "working_no_pod_mpmd" "${LONG_TIMEOUT}" \
        srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
        "${TEST_DIR}/mpmd_working_no_pod.conf"
else
    skip_step "working_no_pod_mpmd" "No working control binary available."
fi

if [ -n "${NODEVICE_NEKO_EXE:-}" ]; then
    run_neko_step "nodevice_no_pod_mpmd" "${LONG_TIMEOUT}" \
        srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
        "${TEST_DIR}/mpmd_nodevice_no_pod.conf"
else
    skip_step "nodevice_no_pod_mpmd" "No no-device-MPI binary available."
fi

if [ "${STEP_RC_MAP[current_no_pod_mpmd]:-1}" = "0" ]; then
    run_neko_step "current_pod_mpmd" "${LONG_TIMEOUT}" \
        srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
        "${TEST_DIR}/mpmd_current_pod.conf"
else
    skip_step "current_pod_mpmd" \
        "Skipped because current_no_pod_mpmd did not pass."
fi

if [ "${STEP_RC_MAP[working_no_pod_mpmd]:-1}" = "0" ]; then
    run_neko_step "working_pod_mpmd" "${LONG_TIMEOUT}" \
        srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
        "${TEST_DIR}/mpmd_working_pod.conf"
else
    skip_step "working_pod_mpmd" \
        "Skipped because working_no_pod_mpmd did not pass."
fi

if [ "${STEP_RC_MAP[nodevice_no_pod_mpmd]:-1}" = "0" ]; then
    run_neko_step "nodevice_pod_mpmd" "${LONG_TIMEOUT}" \
        srun -n "${TOTAL_RANKS}" --label --unbuffered --multi-prog \
        "${TEST_DIR}/mpmd_nodevice_pod.conf"
else
    skip_step "nodevice_pod_mpmd" \
        "Skipped because nodevice_no_pod_mpmd did not pass."
fi

summary_file="${TEST_DIR}/summary.txt"
print_section "SUMMARY"
{
    printf '%-32s %s\n' "step" "rc"
    for idx in "${!STEP_NAMES[@]}"; do
        printf '%-32s %s\n' "${STEP_NAMES[$idx]}" "${STEP_RCS[$idx]}"
    done
} > "${summary_file}"
cat "${summary_file}"

printf "\nArtifacts written to %s\n" "${TEST_DIR}"
printf "This debug example always exits successfully so the framework archives the logs.\n"

exit 0
