#!/usr/bin/env bash
set -e

[ -d "$HOME/miniconda" ] && CONDA_BASE=$HOME/miniconda
[ -d "$HOME/miniconda3" ] && CONDA_BASE=$HOME/miniconda3
[ -f `which conda` ] && CONDA_BASE="$(conda info --base)"

if [ -z "$CONDA_BASE" ]
then
    echo "Could not locate Miniconda in PATH or default directories."
    exit 1
else
    echo "Using Miniconda installation: $CONDA_BASE"
fi

BARCODES_ENV=barcodes
${CONDA_BASE}/bin/conda env create -f env.yml
source "${CONDA_BASE}/bin/activate" ${BARCODES_ENV}
cat > "${HOME}/.barcodesrc" <<EOF
source "${CONDA_BASE}/bin/activate" ${BARCODES_ENV}
EOF
