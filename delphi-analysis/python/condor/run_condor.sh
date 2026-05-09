#!/bin/bash
set -euo pipefail

echo ">>> Starting Condor job at $(date)"
echo ">>> Host: $(hostname)"
echo ">>> Python: $(python3 --version) at $(which python3)"

# --- Move into Condor scratch directory (safe, isolated area)
WORKDIR=${TMPDIR:-$PWD}
cd "$WORKDIR"
echo ">>> Running in scratch: $WORKDIR"

# --- Copy analysis code + this script into scratch
CODEDIR=/afs/cern.ch/user/z/zhangj/private/DELPHI/delphi-nanoaod/delphi-analysis/python
cp $CODEDIR/*.py .
cp $CODEDIR/condor/run_condor.sh .
cp $CODEDIR/efficiency_*.root .
source $CODEDIR/RooUnfold/build/setup.sh

# --- Ensure uproot/awkward dependencies are available ---
python3 -c "import typing_extensions" 2>/dev/null || pip install --user --quiet typing_extensions
python3 -c "import awkward" 2>/dev/null || pip install --user --quiet awkward uproot
export PYTHONPATH=$HOME/.local/lib/python3.9/site-packages:${PYTHONPATH:-}

echo ">>> Files in scratch:"
ls -lh

# --- Show arguments
echo ">>> Arguments:"
echo "  Script: $1"
echo "  Input : $2"
echo "  Output: $3"
echo "  Extra args: ${@:4}"

# --- Run Python
python3 "$1" "$2" "$3" "${@:4}"

 # --- Move result safely to EOS
OUTDIR=/eos/user/z/zhangj/DELPHI/Ntuples/
#OUTDIR=$CODEDIR/out
mkdir -p "$OUTDIR"
mv -v "$3" "$OUTDIR/"

echo ">>> Job finished at $(date)"