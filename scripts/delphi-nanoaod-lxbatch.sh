#!/bin/bash -xe

# parse arguments
if [ "$4" = "default" ]; then
    CONFIG="delphi-nanoaod"
    OUTPUT="$1"
else
    CONFIG="delphi-nanoaod-$4"
    OUTPUT="$1-$4"
fi
if [ "$2" = "." ]; then
    OUTPUT="$OUTPUT.root"
    NICKNAME="$1"
else
    OUTPUT="$OUTPUT-$3.root"
    NICKNAME="$1/$2"
fi
OPTIONS="$5"

# Set up the environment
set +x
if [ -z "$DELPHI" ]; then
    source /cvmfs/delphi.cern.ch/setup.sh
fi
if [ -z "$ROOTSYS" ]; then
    source /cvmfs/sft.cern.ch/lcg/app/releases/ROOT/6.34.08/x86_64-almalinux8.10-gcc85-opt/bin/thisroot.sh
fi
set -x

CONFIG_FILE="$(realpath "$SUBMIT_DIR/../config/$CONFIG.yaml")"
BINARY_FILE="$(realpath "$SUBMIT_DIR/../build/delphi-nanoaod/delphi-nanoaod")"


"$BINARY_FILE" -N "$NICKNAME" -C "$CONFIG_FILE" -O "$OUTPUT" "$OPTIONS"

# ls -l

# xrdcp --force --retry 3 --cksum adler32 "$OUTPUT" "root://eosproject.cern.ch//eos/project/d/delphi/public/nanoaod/v0.2.0/$OUTPUT"
