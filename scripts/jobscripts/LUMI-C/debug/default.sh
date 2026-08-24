#!/bin/bash -l

# LUMI-C submit settings for the split-07 debug example.

# Queue name
#SBATCH --partition=standard

# CPU-only debug layout on LUMI-C.
# Keep the same 8 Neko + 48 Python rank shape used on LUMI-G, but without any
# GPU allocation or GPU-aware MPI.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --cpus-per-task=1

# Time specifications (dd-hh:mm:ss)
#SBATCH --time=00-01:00:00

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

source functions.sh

export OMP_NUM_THREADS="${OMP_NUM_THREADS:-1}"
export MPICH_GPU_SUPPORT_ENABLED=0
export MPICH_DMAPP_APP_IS_WORLD=1
export ATP_ENABLED=true
unset ROCR_VISIBLE_DEVICES

export CASE_FILE="${CASE_FILE:-debug.case}"
export PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python)}"
export PYTHON_SCRIPT="${PYTHON_SCRIPT:-${MAIN_DIR}/scripts/python/pod_state_recover.py}"
export NEKO_RANKS="${NEKO_RANKS:-8}"
export PY_RANKS="${PY_RANKS:-48}"
export DEBUG_USE_GPU=0
export DEBUG_CPUS_PER_TASK="${DEBUG_CPUS_PER_TASK:-1}"
export DEBUG_TIMEOUT_SHORT="${DEBUG_TIMEOUT_SHORT:-120}"
export DEBUG_TIMEOUT_LONG="${DEBUG_TIMEOUT_LONG:-180}"
export DEBUG_TIMEOUT_ATTEMPT="${DEBUG_TIMEOUT_ATTEMPT:-60}"
export DEBUG_TIMEOUT_POD_ATTEMPT="${DEBUG_TIMEOUT_POD_ATTEMPT:-90}"

run "$example"
