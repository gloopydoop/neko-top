#!/usr/bin/env bash
set -euo pipefail

N=10
LX=2
for arg in "$@"; do
    case "${arg}" in
        -N*) N="${arg#-N}" ;;
        -L*) LX="${arg#-L}" ;;
        --Lx=*) LX="${arg#--Lx=}" ;;
        -h|--help)
            echo "Usage: prepare.sh [-N#] [-L#|--Lx=#]"
            exit 0
            ;;
        *) echo "Invalid option: ${arg}" >&2; exit 1 ;;
    esac
done

if [ -z "${NEKO_DIR:-}" ]; then
    MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
    if [ -d "${MAIN_DIR}/external/neko/bin" ]; then
        export NEKO_DIR="${MAIN_DIR}/external/neko"
    fi
fi
if [ -n "${NEKO_DIR:-}" ]; then
    export PATH="${NEKO_DIR}/bin:${PATH}"
    export LD_LIBRARY_PATH="${NEKO_DIR}/lib:${LD_LIBRARY_PATH:-}"
fi
if [ -n "${MAIN_DIR:-}" ]; then
    for dep in neko hdf5 json-fortran adios2; do
        dep_lib="${MAIN_DIR}/external/${dep}/lib"
        if [ -d "${dep_lib}" ]; then
            export LD_LIBRARY_PATH="${dep_lib}:${LD_LIBRARY_PATH:-}"
        fi
    done
fi

command -v genmeshbox >/dev/null 2>&1 || {
    echo "Neko tool 'genmeshbox' not found; set NEKO_DIR or PATH." >&2
    exit 1
}

NX=$((N * LX))
Z=$(python3 -c "print(1.0 / ${N})")
genmeshbox 0 "${LX}" 0 1 0 "${Z}" "${NX}" "${N}" 1 .false. .true. .true.
