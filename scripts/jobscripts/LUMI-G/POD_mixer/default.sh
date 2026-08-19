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

export MPICH_GPU_SUPPORT_ENABLED=1
export ATP_ENABLED=true
export NEKO_STARTUP_DELAY=${NEKO_STARTUP_DELAY:-20}
export NEKO_RANKS_PER_NODE=${NEKO_RANKS_PER_NODE:-8}
export PY_RANKS_PER_NODE=${PY_RANKS_PER_NODE:-48}
export SRUN_MPMD_FLAGS="${SRUN_MPMD_FLAGS:---distribution=block:block}"
export POD_MIXER_MPMD_CONF_GENERATOR="${MAIN_DIR}/scripts/jobscripts/LUMI-G/POD_mixer/generate_mpmd_conf.sh"

run $example
