#!/bin/bash

# Read positional arguments
INPUT_FILE="$1"
OUTPUT_FILE="$2"
TPCNTUPLE_DIR_PATH="$3"
DELPHINANOAOD_DIR_PATH="$4"
DATA_OR_MC="$5"

# Print the arguments
echo "Arguments passed to this script:"
echo "  input file         : $INPUT_FILE"
echo "  output file        : $OUTPUT_FILE"
echo "  TPCNtuple dir      : $TPCNTUPLE_DIR_PATH"
echo "  DelphiNanoAOD dir  : $DELPHINANOAOD_DIR_PATH"
echo "  data or mc         : $DATA_OR_MC"

# Create a temporary working directory
WORKDIR=$(mktemp -d /tmp/delphi_job_XXXXXX)
echo "Working directory: $WORKDIR"
trap 'rm -rf "$WORKDIR"' EXIT

# Copy everything into the temp directory
BASEDIR="/afs/cern.ch/user/z/zhangj/private/DELPHI/delphi-nanoaod"
cp "$BASEDIR/build/delphi-nanoaod/delphi-nanoaod" "$WORKDIR/"
cp "$BASEDIR/config/delphi-nanoaod.yaml" "$WORKDIR/"
cp "$BASEDIR/scripts/treefy.C" "$WORKDIR/"
cd "$WORKDIR" || exit 1

# Set up environments
source /cvmfs/delphi.cern.ch/setup.sh
source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.34.04/x86_64-almalinux9.5-gcc115-opt/bin/thisroot.sh

# Write unique PDLINPUT file
echo "FILE=$INPUT_FILE" > "${OUTPUT_FILE}_dummy"

# Run nanoAOD producer
if [ "$DATA_OR_MC" = "MC" ]; then
    delphi-nanoaod -P "${OUTPUT_FILE}_dummy" --mc --config delphi-nanoaod.yaml --output "$OUTPUT_FILE"
else
    delphi-nanoaod -P "${OUTPUT_FILE}_dummy" --config delphi-nanoaod.yaml --output "$OUTPUT_FILE"
fi

# Run treefy step
root -q -b -l "treefy.C+(\"$OUTPUT_FILE\")"

# Route output files based on their names
tpc="$OUTPUT_FILE"
nanotree="${tpc/.root/_ttree.root}"

# Move files to appropriate directories
echo "Moving TPCNtuple file to: $TPCNTUPLE_DIR_PATH"
mv "$tpc" "$TPCNTUPLE_DIR_PATH/"

echo "Moving DelphiNanoAOD file to: $DELPHINANOAOD_DIR_PATH"
mv "$nanotree" "$DELPHINANOAOD_DIR_PATH/"

echo "Job completed successfully!"