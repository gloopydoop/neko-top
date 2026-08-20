#!/bin/bash -l

# LUMI-G submit settings for the POD mixer example. This case uses a full-node
# MPMD layout: eight GPU-backed Neko ranks plus forty-eight CPU-only Python
# ranks per node.

# Queue name
#SBATCH --partition=standard-g

# Full-node POD layout on LUMI-G.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --gpus-per-node=8
#SBATCH --mem=480GB

# Time specifications (dd-hh:mm:ss)
#SBATCH --time=00-00:30:00

# Notification options
#SBATCH --mail-type=END

# Mandatory options
#SBATCH --open-mode=append
#SBATCH --output=output.log
#SBATCH --error=error.log

set -e

if [[ -z "$SLURM_JOB_NAME" && (($# > 0)) ]]; then
    example=$1
elif [ -n "${SLURM_JOB_NAME:-}" ]; then
    example=$SLURM_JOB_NAME
else
    printf "ERROR: No example supplied" >&2
    exit 1
fi

ml craype-accel-amd-gfx90a rocm

cat <<'EOF' > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
exec "$@"
EOF

chmod +x ./select_gpu

source functions.sh
source "${MAIN_DIR}/scripts/mpmd_run_helpers.sh"

export MPICH_GPU_SUPPORT_ENABLED=1
# Cray MPICH requires this for srun --multi-prog MPMD jobs on LUMI-G.
export MPICH_DMAPP_APP_IS_WORLD=1
export ATP_ENABLED=true
export NEKO_STARTUP_DELAY=${NEKO_STARTUP_DELAY:-20}
export NEKO_RANKS_PER_NODE=${NEKO_RANKS_PER_NODE:-8}
export PY_RANKS_PER_NODE=${PY_RANKS_PER_NODE:-48}
export POD_MIXER_DEBUG=${POD_MIXER_DEBUG:-1}
trap 'rm -f ./select_gpu ./mpmd.conf' EXIT

function pod_mixer_write_debug_snapshot() {
    local casefile=$1
    local logfile=$2
    local neko_exe=$3
    local py_script=$4
    local debug_log=./mpmd_runtime_debug.log
    local rank_probe_log=./mpmd_rank_probe.log
    local select_gpu_probe_log=./mpmd_select_gpu_probe.log

    cp -f ./mpmd.conf ./mpmd.conf.debug
    cp -f ./select_gpu ./select_gpu.debug.sh

    {
        printf "Timestamp: %s\n" "$(date -Is)"
        printf "Hostname: %s\n" "$(hostname)"
        printf "PWD: %s\n" "$(pwd)"
        printf "Case file: %s\n" "${casefile}"
        printf "Case path: %s\n" "$(realpath "${casefile}")"
        printf "Log file: %s\n" "${logfile}"
        printf "Neko exe: %s\n" "${neko_exe}"
        printf "Neko exe path: %s\n" "$(realpath "${neko_exe}")"
        printf "Python bin: %s\n" "${PYTHON_BIN}"
        printf "Python script: %s\n" "${py_script}"
        printf "Python script path: %s\n" "$(realpath "${py_script}")"
        printf "srun path: %s\n" "$(command -v srun 2>/dev/null || echo '<not found>')"
        printf "mpicc path: %s\n" "$(command -v mpicc 2>/dev/null || echo '<not found>')"
        printf "sbatch path: %s\n" "$(command -v sbatch 2>/dev/null || echo '<not found>')"
        if command -v srun >/dev/null 2>&1; then
            printf "srun version: %s\n" "$(srun --version 2>&1 | head -n 1)"
        fi
        printf -- "---- module list ----\n"
        if command -v module >/dev/null 2>&1; then
            module list 2>&1 || true
        else
            printf "module command not available\n"
        fi
        printf -- "---- selected environment ----\n"
        env | sort | grep -E \
            '^(SLURM|PMI|PMIX|MPICH|ROCR|HIP|CUDA|LD_LIBRARY_PATH|PYTHONPATH|PATH|VIRTUAL_ENV|CONDA_PREFIX|ADIOS2|NEKO|FI_|UCX_|OMP_)=' || true
        printf -- "---- scontrol show job ----\n"
        if [ -n "${SLURM_JOB_ID:-}" ] && command -v scontrol >/dev/null 2>&1; then
            scontrol show job -d "${SLURM_JOB_ID}" 2>&1 || true
        else
            printf "scontrol or SLURM_JOB_ID unavailable\n"
        fi
        printf -- "---- ldd neko ----\n"
        ldd "${neko_exe}" 2>&1 || true
        printf -- "---- python probe ----\n"
        "${PYTHON_BIN}" - <<'EOF'
import os
import sys

print(f"sys.executable: {sys.executable}")
print(f"sys.version: {sys.version.replace(chr(10), ' ')}")
for name in ("mpi4py", "adios2"):
    try:
        mod = __import__(name)
        print(f"{name} file: {getattr(mod, '__file__', '<builtin>')}")
        print(f"{name} version: {getattr(mod, '__version__', '<unknown>')}")
    except Exception as exc:
        print(f"{name} import failed: {exc!r}")

try:
    from mpi4py import MPI
    print(f"MPI vendor: {MPI.get_vendor()}")
    print("MPI library version:")
    print(MPI.Get_library_version().strip())
except Exception as exc:
    print(f"mpi4py.MPI probe failed: {exc!r}")

for key in (
    "LD_LIBRARY_PATH",
    "PYTHONPATH",
    "VIRTUAL_ENV",
    "CONDA_PREFIX",
    "NEKO_COMM_ID",
    "NEKO_CTRL_PEER_ROOT",
):
    print(f"{key}: {os.getenv(key, '<unset>')}")
EOF
        printf -- "---- mpmd.conf ----\n"
        cat ./mpmd.conf
        printf -- "---- select_gpu ----\n"
        cat ./select_gpu
    } > "${debug_log}" 2>&1

    if command -v srun >/dev/null 2>&1; then
        {
            srun --unbuffered -n "${SLURM_NTASKS:-56}" /bin/bash -lc \
                'printf "rank=%s local=%s node=%s host=%s pmi_rank=%s pmi_size=%s pmix_rank=%s pmix_size=%s rocr=%s\n" \
                "${SLURM_PROCID:-<unset>}" "${SLURM_LOCALID:-<unset>}" \
                "${SLURM_NODEID:-<unset>}" "$(hostname)" "${PMI_RANK:-<unset>}" \
                "${PMI_SIZE:-<unset>}" "${PMIX_RANK:-<unset>}" \
                "${PMIX_SIZE:-<unset>}" "${ROCR_VISIBLE_DEVICES:-<unset>}"'
        } 2>&1 | sort -V > "${rank_probe_log}" || true

        {
            srun --unbuffered -n "${NEKO_RANKS_PER_NODE}" ./select_gpu /bin/bash -lc \
                'printf "rank=%s local=%s host=%s rocr=%s pmi_rank=%s pmi_size=%s\n" \
                "${SLURM_PROCID:-<unset>}" "${SLURM_LOCALID:-<unset>}" \
                "$(hostname)" "${ROCR_VISIBLE_DEVICES:-<unset>}" \
                "${PMI_RANK:-<unset>}" "${PMI_SIZE:-<unset>}"'
        } 2>&1 | sort -V > "${select_gpu_probe_log}" || true
    fi

    printf "Wrote debug files:\n"
    printf "\t%s\n" "${debug_log}"
    printf "\t%s\n" "${rank_probe_log}"
    printf "\t%s\n" "${select_gpu_probe_log}"
    printf "\t%s\n" "./mpmd.conf.debug"
    printf "\t%s\n" "./select_gpu.debug.sh"
}

function run {
    set +e

    casefile=($(find . -name "*.case" -or -name "*.json"))
    if [[ ${#casefile[@]} -eq 0 ]]; then
        printf >&2 "ERROR: No case file found.\n"
        return 1
    elif [[ ${#casefile[@]} -eq 1 ]]; then
        casefile=${casefile[0]}
        logfile=$(basename -- ${casefile%.*}).log
    else
        logfile=$(basename -- $(dirname $(realpath $0))).log
    fi

    if [ -s error.log ]; then
        grep "ERROR: Optimizer stopped after reaching the maximum runtime" \
            error.log >/dev/null || return 1
        echo "" > error.log
    fi

    if [ -f "$logfile" ]; then
        old_run=run_$(find ./ -maxdepth 1 -type d -name "run_*" | wc -l)
        old_run=$(printf "%s_%02d" "run" $((10#${old_run#run_} + 1)))
        mkdir -p ./$old_run

        find ./ -maxdepth 1 -not -empty -type f -name "*.log" \
            -not -name "output.log" -not -name "error.log" \
            -exec mv -ft ./$old_run {} \;
        cp -ft ./$old_run output.log
        [ -s error.log ] && cp -ft ./$old_run error.log

        printf "Ready" >./output.log
    fi

    printf "Executing Neko.\n" > ./output.log
    printf "See $logfile for the status output.\n"
    export NEKO_LOG_FILE=$logfile

    prepare 2>error.log || return 1
    if [ -s ./error.log ]; then
        printf "ERROR: An error occured during preparation.\n"
        printf "See error.log for details.\n"
        return 1
    fi
    rm -fr error.log && touch error.log

    printf "=%.0s" {1..80} && printf "\n"
    printf "Running example: %s.\n" $example
    printf "=%.0s" {1..80} && printf "\n"

    generator="${MAIN_DIR}/scripts/jobscripts/LUMI-G/POD_mixer/generate_mpmd_conf.sh"
    if [ ! -x "${generator}" ]; then
        echo "Error: missing POD_mixer MPMD generator: ${generator}" >&2
        return 1
    fi

    py_script="${PYTHON_SCRIPT:-${MAIN_DIR}/scripts/python/pod_state_recover.py}"
    neko_exe="${neko}"
    if [[ "${neko_exe}" == "./select_gpu "* ]]; then
        neko_exe="${neko_exe#./select_gpu }"
    fi
    mpmd_prepare_python_runtime "${MAIN_DIR}" || return 1
    mpmd_print_runtime_env

    "${generator}" \
        --nodes "${SLURM_NNODES:-1}" \
        --case-file "${casefile}" \
        --neko-exe "${neko_exe}" \
        --python-bin "${PYTHON_BIN}" \
        --python-script "${py_script}" \
        --neko-ranks-per-node "${NEKO_RANKS_PER_NODE}" \
        --python-ranks-per-node "${PY_RANKS_PER_NODE}" \
        --output "./mpmd.conf" \
        --select-gpu "./select_gpu" || return 1

    if [ "${POD_MIXER_DEBUG}" != "0" ]; then
        pod_mixer_write_debug_snapshot "${casefile}" "${logfile}" \
            "${neko_exe}" "${py_script}"
    fi

    TIME_START=$(date +%s)
    if [ -n "${SRUN_MPMD_FLAGS:-}" ]; then
        read -r -a srun_extra <<< "${SRUN_MPMD_FLAGS}"
        srun --unbuffered "${srun_extra[@]}" --multi-prog ./mpmd.conf \
            > "${logfile}" 2>error.log
    else
        srun --unbuffered --multi-prog ./mpmd.conf > "${logfile}" 2>error.log
    fi
    rc=$?
    TIME_END=$(date +%s)

    if [ "${rc}" -ne 0 ]; then
        printf "ERROR: An error occurred during execution.\n"
        printf "See error.log for details.\n"
        return 1
    fi

    if [ -s ./error.log ]; then
        printf "ERROR: An error occurred during execution.\n"
        printf "See error.log for details.\n"
        return 1
    fi

    printf "\nExample concluded.\n"
    TIME_DIFF=$((TIME_END - TIME_START))
    printf "Execution time: %02d:%02d:%02d\n" \
        $((TIME_DIFF / 3600)) $((TIME_DIFF % 3600 / 60)) $((TIME_DIFF % 60))

    cleanup 2>error.log || return 1

    if [ -s ./error.log ]; then
        printf "ERROR: An error occurred during cleanup.\n"
        printf "See error.log for details.\n"
        return 1
    fi
}

run $example
