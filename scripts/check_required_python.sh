#!/usr/bin/env bash

# 1. Check Python version
if ! command -v python3 &>/dev/null; then
    echo "ERROR: Python3 is not installed." >&2
    exit 1
fi

PY_VER=$(python3 -c 'import sys; print(f"{sys.version_info.major}.{sys.version_info.minor}")')
REQ_VER="3.8"
if [[ "$(printf '%s\n' "$REQ_VER" "$PY_VER" | sort -V | head -n1)" != "$REQ_VER" ]]; then
    echo "ERROR: Python $REQ_VER or higher is required. Found $PY_VER." >&2
    exit 1
fi

# 2. Check required modules
modules=(
    numpy
    pysam
    pandas
    matplotlib
    h5py
    scipy
)

missing=()
for mod in "${modules[@]}"; do
    python3 -c "import $mod" 2>/dev/null || missing+=("$mod")
done

if (( ${#missing[@]} )); then
    echo "ERROR: Missing Python modules:" >&2
    for m in "${missing[@]}"; do
        echo "  - $m" >&2
    done
    echo -e "\nInstall them with:" >&2
    echo "  pip install ${missing[*]}" >&2
    exit 1
fi

echo "Python version: $PY_VER"
echo "Required modules are installed: ${modules[*]}"
exit 0
