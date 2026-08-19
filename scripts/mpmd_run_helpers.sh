#!/bin/bash

# Generic helpers for coupled Neko/Python MPMD launches.
#
# This layer only covers:
# - locating the Python runtime recorded by setup.sh
# - validating that runtime before launching
# - simple mixed launches via mpirun or srun --multi-prog
#
# Cluster-specific placement policies are intentionally kept out of this file.
# If scripts/mpmd_slurm_helpers.sh exists, mpmd_launch_shared() will delegate
# to its mpmd_launch_shared_bound() entry point when that helper declares
# itself applicable.

_mpmd_helper_dir=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)
if [ -f "${_mpmd_helper_dir}/mpmd_slurm_helpers.sh" ]; then
    # shellcheck disable=SC1091
    source "${_mpmd_helper_dir}/mpmd_slurm_helpers.sh"
fi
unset _mpmd_helper_dir

function mpmd_python_exe() {
    if [ -n "${PYTHON_BIN:-}" ]; then
        if [ -x "${PYTHON_BIN}" ]; then
            printf '%s\n' "${PYTHON_BIN}"
            return 0
        fi

        if command -v "${PYTHON_BIN}" >/dev/null 2>&1; then
            command -v "${PYTHON_BIN}"
            return 0
        fi

        echo "Error: PYTHON_BIN is set but not executable: ${PYTHON_BIN}" >&2
        return 1
    fi

    command -v python3 2>/dev/null || command -v python 2>/dev/null
}

function mpmd_find_repo_root() {
    local start_dir=$1
    local root_dir

    root_dir=$(realpath "${start_dir}")

    while [ ! -d "${root_dir}/sources" ] && [ "${root_dir}" != "/" ]; do
        root_dir="$(dirname "${root_dir}")"
    done

    if [ ! -d "${root_dir}/sources" ]; then
        echo "Error: could not locate repo root from ${start_dir}" >&2
        return 1
    fi

    printf '%s\n' "${root_dir}"
}

function mpmd_runtime_env_file() {
    local root_dir=$1
    local repo_root

    repo_root=$(mpmd_find_repo_root "${root_dir}") || return 1
    printf '%s\n' "${repo_root}/build/mpmd_runtime.env"
}

function mpmd_source_runtime_env() {
    local root_dir=$1
    local runtime_env

    runtime_env=$(mpmd_runtime_env_file "${root_dir}") || return 1
    if [ ! -f "${runtime_env}" ]; then
        echo "Error: Python runtime env not found: ${runtime_env}" >&2
        echo "Run ./setup.sh -e from the repo root after activating the" >&2
        echo "target Python environment." >&2
        return 1
    fi

    # shellcheck disable=SC1090
    source "${runtime_env}"
    export MPMD_RUNTIME_ENV_FILE="${runtime_env}"
}

function mpmd_validate_python_runtime() {
    local root_dir=$1
    local repo_root
    local validator
    local pyexe

    repo_root=$(mpmd_find_repo_root "${root_dir}") || return 1
    validator="${repo_root}/scripts/python/validate_mpmd_runtime.py"
    if [ ! -f "${validator}" ]; then
        echo "Error: Python runtime validator not found: ${validator}" >&2
        return 1
    fi

    pyexe=$(mpmd_python_exe) || {
        echo "Error: could not find python3 or python in PATH." >&2
        return 1
    }

    if ! "${pyexe}" "${validator}"; then
        echo "Error: the active Python runtime is not valid." >&2
        echo "Python: ${pyexe}" >&2
        echo "Runtime env: ${MPMD_RUNTIME_ENV_FILE:-<unset>}" >&2
        return 1
    fi

    export PYTHON_BIN="${pyexe}"
}

function mpmd_prepare_python_runtime() {
    local root_dir=$1

    mpmd_source_runtime_env "${root_dir}" || return 1
    mpmd_validate_python_runtime "${root_dir}" || return 1
}

function mpmd_ensure_python_runtime() {
    local root_dir=${1:-${MAIN_DIR:-.}}
    mpmd_prepare_python_runtime "${root_dir}"
}

function mpmd_ensure_adios2_python() {
    mpmd_ensure_python_runtime "$@"
}

function mpmd_startup_delay() {
    if [ -n "${NEKO_STARTUP_DELAY:-}" ]; then
        printf '%s\n' "${NEKO_STARTUP_DELAY}"
        return 0
    fi

    printf '%s\n' "${MPMD_NEKO_STARTUP_DELAY:-20}"
}

function mpmd_selected_launcher() {
    local requested=${NEKO_LAUNCHER:-auto}

    case "${requested}" in
        auto)
            if [ -n "${SLURM_JOB_ID:-}" ] && command -v srun >/dev/null 2>&1
            then
                printf 'srun\n'
                return 0
            fi

            if command -v mpirun >/dev/null 2>&1; then
                printf 'mpirun\n'
                return 0
            fi

            if command -v srun >/dev/null 2>&1; then
                printf 'srun\n'
                return 0
            fi

            echo "Error: neither mpirun nor srun was found in PATH." >&2
            return 1
            ;;
        srun|mpirun)
            if ! command -v "${requested}" >/dev/null 2>&1; then
                echo "Error: ${requested} was requested via" >&2
                echo "NEKO_LAUNCHER but is not available in PATH." >&2
                return 1
            fi

            printf '%s\n' "${requested}"
            return 0
            ;;
        *)
            echo "Error: NEKO_LAUNCHER must be one of:" >&2
            echo "auto, srun, mpirun." >&2
            return 1
            ;;
    esac
}

function mpmd_print_runtime_env() {
    echo "Using Python:       $(mpmd_python_exe 2>/dev/null || echo \
'<not found>')"
    echo "Using mpicc:        $(command -v mpicc 2>/dev/null || echo \
'<not found>')"
    echo "Using mpirun:       $(command -v mpirun 2>/dev/null || echo \
'<not found>')"
    echo "Using srun:         $(command -v srun 2>/dev/null || echo \
'<not found>')"
    echo "NEKO_LAUNCHER:      ${NEKO_LAUNCHER:-auto}"
    echo "Selected launcher:  $(mpmd_selected_launcher 2>/dev/null || echo \
'<unavailable>')"
    echo "Python runtime env: ${MPMD_RUNTIME_ENV_FILE:-<unset>}"
    echo "CONDA_PREFIX:       ${CONDA_PREFIX:-<unset>}"
    echo "VIRTUAL_ENV:        ${VIRTUAL_ENV:-<unset>}"
    echo "ADIOS2_PATH:        ${ADIOS2_PATH:-<unset>}"
    echo "PYTHONPATH:         ${PYTHONPATH:-<unset>}"
    echo "LD_LIBRARY_PATH:    ${LD_LIBRARY_PATH:-<unset>}"
    echo "---------------------------------------------------------"
}

function mpmd_resolve_case_file() {
    local script_dir=$1
    local requested=${2:-}
    local case_files

    if [ -n "${requested}" ]; then
        printf '%s\n' "${requested}"
        return 0
    fi

    shopt -s nullglob
    case_files=( "${script_dir}"/*.case )
    shopt -u nullglob

    if [ ${#case_files[@]} -eq 1 ]; then
        printf '%s\n' "$(basename "${case_files[0]}")"
        return 0
    fi

    if [ ${#case_files[@]} -eq 0 ]; then
        echo "Error: no .case file found in ${script_dir}" >&2
    else
        echo "Error: multiple .case files found in ${script_dir}" >&2
    fi
    return 1
}

function mpmd_launch_shared() {
    local case_path=$1
    local py_script=$2
    local neko_exe=$3
    local py_ranks=$4
    local neko_ranks=$5
    local log_file=$6

    local total_ranks=$((py_ranks + neko_ranks))
    local conf_file
    local abs_case_path
    local abs_case_dir
    local abs_py_script
    local abs_neko_exe
    local runtime_env
    local pyexe
    local startup_delay
    local launcher
    local errexit_was_set
    local rc
    local i
    local neko_cmd
    local py_cmd
    local -a mpirun_cmd

    abs_case_path=$(realpath "${case_path}")
    abs_case_dir=$(dirname "${abs_case_path}")
    abs_py_script=$(realpath "${py_script}")
    abs_neko_exe=$(realpath "${neko_exe}")
    runtime_env=${MPMD_RUNTIME_ENV_FILE:-}
    if [ -z "${runtime_env}" ]; then
        runtime_env=$(mpmd_runtime_env_file "${abs_case_dir}") || return 1
    fi
    runtime_env=$(realpath "${runtime_env}")
    launcher=$(mpmd_selected_launcher) || return 1

    if declare -F mpmd_launch_shared_bound_available >/dev/null 2>&1 &&
        mpmd_launch_shared_bound_available
    then
        mpmd_launch_shared_bound \
            "${case_path}" "${py_script}" "${neko_exe}" \
            "${py_ranks}" "${neko_ranks}" "${log_file}"
        return $?
    fi

    pyexe=$(mpmd_python_exe) || {
        echo "Error: could not find python in PATH." >&2
        return 1
    }
    startup_delay=$(mpmd_startup_delay)

    if [ "${launcher}" = "srun" ]; then
        conf_file=$(mktemp "${TMPDIR:-/tmp}/mpmd.XXXXXX.conf")

        {
            for ((i=0; i<py_ranks; i++)); do
                echo "${i} /bin/bash -c 'cd \"${abs_case_dir}\" &&" \
                    "source \"${runtime_env}\" &&" \
                    "exec /usr/bin/env NEKO_COMM_ID=1" \
                    "NEKO_CTRL_PEER_ROOT=${py_ranks} \"${pyexe}\"" \
                    "\"${abs_py_script}\" \"${abs_case_path}\"'"
            done

            for ((i=py_ranks; i<total_ranks; i++)); do
                echo "${i} /bin/bash -c 'cd \"${abs_case_dir}\" &&" \
                    "source \"${runtime_env}\" &&" \
                    "sleep ${startup_delay};" \
                    "exec /usr/bin/env NEKO_COMM_ID=0" \
                    "NEKO_CTRL_PEER_ROOT=0 \"${abs_neko_exe}\"" \
                    "\"${abs_case_path}\"'"
            done
        } > "${conf_file}"

        echo "Running MPMD job with ${total_ranks} ranks"
        echo "Launcher: srun --multi-prog"
        echo "Config:   ${conf_file}"
        echo "Output:   ${log_file}"

        errexit_was_set=0
        if [[ $- == *e* ]]; then
            errexit_was_set=1
            set +e
        fi

        srun --unbuffered --multi-prog "${conf_file}" \
            > "${log_file}" 2>&1
        rc=$?

        if [ "${errexit_was_set}" -eq 1 ]; then
            set -e
        fi

        if [ "${rc}" -ne 0 ]; then
            echo "Error: shared MPMD launch failed. See ${log_file}." >&2
        fi

        rm -f "${conf_file}"
        return ${rc}
    fi

    printf -v py_cmd \
        'cd %q && source %q && exec /usr/bin/env NEKO_COMM_ID=1 ' \
        "${abs_case_dir}" "${runtime_env}"
    printf -v py_cmd \
        '%sNEKO_CTRL_PEER_ROOT=%q %q %q %q' \
        "${py_cmd}" "${py_ranks}" "${pyexe}" "${abs_py_script}" \
        "${abs_case_path}"

    printf -v neko_cmd \
        'cd %q && source %q && sleep %q; exec /usr/bin/env NEKO_COMM_ID=0 ' \
        "${abs_case_dir}" "${runtime_env}" "${startup_delay}"
    printf -v neko_cmd \
        '%sNEKO_CTRL_PEER_ROOT=0 %q %q' \
        "${neko_cmd}" "${abs_neko_exe}" "${abs_case_path}"

    mpirun_cmd=(
        mpirun
        --tag-output
        -n "${py_ranks}"
        /bin/bash
        -c
        "${py_cmd}"
        :
        -n "${neko_ranks}"
        /bin/bash
        -c
        "${neko_cmd}"
    )

    echo "Launching shared MPI job:"
    printf '  %q' "${mpirun_cmd[@]}"
    printf '\n'
    echo "Output:            ${log_file}"

    errexit_was_set=0
    if [[ $- == *e* ]]; then
        errexit_was_set=1
        set +e
    fi

    "${mpirun_cmd[@]}" > "${log_file}" 2>&1
    rc=$?

    if [ "${errexit_was_set}" -eq 1 ]; then
        set -e
    fi

    if [ "${rc}" -ne 0 ]; then
        echo "Error: shared MPMD launch failed. See ${log_file}." >&2
    fi

    return ${rc}
}
