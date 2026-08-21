#!/bin/bash -l

# LUMI-G submit settings for the split-07 debug example.

# Queue name
#SBATCH --partition=standard-g

# Full-node MPMD debug layout on LUMI-G.
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=56
#SBATCH --gpus-per-node=8
#SBATCH --mem=480GB

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

ml craype-accel-amd-gfx90a rocm

cat <<'EOF' > select_gpu
#!/bin/bash

export ROCR_VISIBLE_DEVICES=${SLURM_LOCALID:-0}
exec "$@"
EOF

chmod +x ./select_gpu
trap 'rm -f ./select_gpu' EXIT

source functions.sh

export MPICH_GPU_SUPPORT_ENABLED=1
export MPICH_DMAPP_APP_IS_WORLD=1
export ATP_ENABLED=true

export CASE_FILE="${CASE_FILE:-debug.case}"
export PYTHON_BIN="${PYTHON_BIN:-$(command -v python3 || command -v python)}"
export PYTHON_SCRIPT="${PYTHON_SCRIPT:-${MAIN_DIR}/scripts/python/pod_state_recover.py}"
export NEKO_RANKS="${NEKO_RANKS:-8}"
export PY_RANKS="${PY_RANKS:-48}"
export NEKO_STARTUP_DELAY="${NEKO_STARTUP_DELAY:-20}"
export DEBUG_TIMEOUT_SHORT="${DEBUG_TIMEOUT_SHORT:-120}"
export DEBUG_TIMEOUT_LONG="${DEBUG_TIMEOUT_LONG:-180}"
export DEBUG_TIMEOUT_ATTEMPT="${DEBUG_TIMEOUT_ATTEMPT:-60}"
export DEBUG_TIMEOUT_POD_ATTEMPT="${DEBUG_TIMEOUT_POD_ATTEMPT:-90}"

run "$example"
