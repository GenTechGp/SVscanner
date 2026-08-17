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

# pysam bundles its own htslib and <0.24.0 resolves INFO/END against INFO/SVLEN
# differently. src/extract_sv.py no longer depends on that (see get_sv_end), so
# this is a warning rather than an error - but the version belongs in the log:
# an unexplained change in extracted coordinates is almost always this.
PYSAM_MIN="0.24.0"

missing=()
declare -A mod_version mod_path

for mod in "${modules[@]}"; do
    if ! info=$(python3 -c "
import importlib
m = importlib.import_module('${mod}')
print(getattr(m, '__version__', 'unknown'))
print(getattr(m, '__file__', 'unknown'))
" 2>/dev/null); then
        missing+=("$mod")
        continue
    fi
    mod_version["$mod"]="${info%%$'\n'*}"
    mod_path["$mod"]="${info#*$'\n'}"
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

echo "Python version: $PY_VER ($(command -v python3))"
echo "Required modules are installed:"
for mod in "${modules[@]}"; do
    printf '  %-12s %-10s %s\n' "$mod" "${mod_version[$mod]}" "${mod_path[$mod]}"
done

if [[ "$(printf '%s\n' "$PYSAM_MIN" "${mod_version[pysam]}" | sort -V | head -n1)" != "$PYSAM_MIN" ]]; then
    echo "Warning: pysam ${mod_version[pysam]} is older than ${PYSAM_MIN}, which SVscanner is tested and deployed with." >&2
    echo "Warning: SV spans are resolved by src/extract_sv.py, so results should be unaffected." >&2
fi

exit 0
