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

source functions.sh
source "${MAIN_DIR}/scripts/mpmd_run_helpers.sh"

export MPICH_GPU_SUPPORT_ENABLED=1
export ATP_ENABLED=true
export NEKO_RANKS_PER_NODE=${NEKO_RANKS_PER_NODE:-8}
export PY_RANKS_PER_NODE=${PY_RANKS_PER_NODE:-48}
trap 'rm -f ./select_gpu ./mpmd.conf' EXIT

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
    mpmd_prepare_python_runtime "${MAIN_DIR}" || return 1
    mpmd_print_runtime_env

    "${generator}" \
        --nodes "${SLURM_NNODES:-1}" \
        --case-file "${casefile}" \
        --neko-exe "${neko}" \
        --python-bin "${PYTHON_BIN}" \
        --python-script "${py_script}" \
        --neko-ranks-per-node "${NEKO_RANKS_PER_NODE}" \
        --python-ranks-per-node "${PY_RANKS_PER_NODE}" \
        --output "./mpmd.conf" \
        --select-gpu "./select_gpu" || return 1

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
