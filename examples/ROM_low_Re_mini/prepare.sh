#!/usr/bin/env bash
set -euo pipefail

N=${N:-4}
LENGTH=${LENGTH:-6}
for arg in "$@"; do
    case "${arg}" in
        -N*) N="${arg#-N}" ;;
        -L*) LENGTH="${arg#-L}" ;;
        --length=*) LENGTH="${arg#--length=}" ;;
        -h|--help)
            echo "Usage: prepare.sh [-N#] [-L#|--length=#]"
            exit 0
            ;;
        *) echo "Invalid option: ${arg}" >&2; exit 1 ;;
    esac
done

if [ -z "${NEKO_DIR:-}" ]; then
    MINI_MAIN_DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")/../.." && pwd)
    if [ -d "${MINI_MAIN_DIR}/external/neko/bin" ]; then
        export NEKO_DIR="${MINI_MAIN_DIR}/external/neko"
    fi
fi
if [ -n "${NEKO_DIR:-}" ]; then
    export PATH="${NEKO_DIR}/bin:${PATH}"
fi

command -v genmeshbox >/dev/null 2>&1 || {
    echo "Neko tool genmeshbox not found; set NEKO_DIR or source prepare.env." >&2
    exit 1
}

NX=$((N * LENGTH))
genmeshbox 0 "${LENGTH}" 0 2 0 2 "${NX}" "${N}" "${N}" .false. .false. .false.
