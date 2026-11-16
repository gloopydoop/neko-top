#!/usr/bin/env bash
set -euo pipefail

N=${N:-6}
LX=${LX:-2}
NX=$((N * LX))
NY=$N
NZ=1
Z=$(python3 -c "print(1.0 / ${N})")

if ! command -v genmeshbox >/dev/null 2>&1; then
    echo "Neko tool 'genmeshbox' not found." >&2
    exit 1
fi

genmeshbox 0 "${LX}" 0 1 0 "${Z}" "${NX}" "${NY}" "${NZ}" \
    .false. .true. .true.
