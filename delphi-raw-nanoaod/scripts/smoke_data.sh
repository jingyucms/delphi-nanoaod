#!/bin/bash
# Data smoke test: feed a DELPHI long-DST (.al) file directly into
# delphi-raw-nanoaod and produce a brief sanity report.
#
# Sibling to end_to_end_smoke.sh: that script exercises the full
# Pythia8 -> DELSIM -> (.sdst + .fadana) -> reader chain for MC and
# needs a cmssw-el9 singularity image + the Delphi-Sim-Pipeline repo.
# For real data none of that is required: the .al long-DST already
# carries both MVDH (VD hits) and PA.TE* (per-detector track elements)
# in the same file, so a single reader pass covers everything.
#
# Runs on the host (lxplus EL9): sources the repo's setup.sh (CVMFS
# DELPHI + CVMFS ROOT 6.34.08), writes a one-line PDL input that
# points at the chosen .al file, runs the reader, then opens the
# resulting RNTuple and prints per-collection sums.
#
# Usage
# -----
#   smoke_data.sh [--al PATH] [--events N] [--out DIR] [--raw-bin PATH]
#
# Defaults
# --------
#   --al       /eos/opendata/delphi/collision-data/Y13718/Y13718.100.al
#   --events   100
#   --out      /tmp/raw_smoke_data
#   --raw-bin  $(this repo)/build/delphi-raw-nanoaod/delphi-raw-nanoaod
#
# Any non-default value can also come from an env var of the same
# uppercased name, e.g. `EVENTS=500 ./smoke_data.sh`.

set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

AL=${AL:-/eos/opendata/delphi/collision-data/Y13718/Y13718.100.al}
EVENTS=${EVENTS:-100}
OUT=${OUT:-/tmp/raw_smoke_data}
RAW_BIN=${RAW_BIN:-$REPO_ROOT/build/delphi-raw-nanoaod/delphi-raw-nanoaod}
JOB_ID=${JOB_ID:-smoke_$(date +%Y%m%d_%H%M%S)}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --al)       AL=$2;      shift 2;;
        --events)   EVENTS=$2;  shift 2;;
        --out)      OUT=$2;     shift 2;;
        --raw-bin)  RAW_BIN=$2; shift 2;;
        --job-id)   JOB_ID=$2;  shift 2;;
        -h|--help)  sed -n '2,30p' "$0"; exit 0;;
        *) echo "Unknown flag: $1" >&2; exit 2;;
    esac
done

echo "=== config ==="
echo "  al       = $AL"
echo "  events   = $EVENTS"
echo "  out      = $OUT"
echo "  raw-bin  = $RAW_BIN"
echo "  job-id   = $JOB_ID"
echo

[[ -f "$AL" ]]      || { echo "missing long-DST file: $AL" >&2; exit 1; }
[[ -x "$RAW_BIN" ]] || {
    echo "missing raw-reader binary: $RAW_BIN" >&2
    echo "build it first, e.g.:" >&2
    echo "  source $REPO_ROOT/setup.sh" >&2
    echo "  cmake -B $REPO_ROOT/build -DROOT_DIR=\"\$ROOTSYS/cmake\" -DCMAKE_PREFIX_PATH=\"\$ROOTSYS\"" >&2
    echo "  cmake --build $REPO_ROOT/build --target delphi-raw-nanoaod -j4" >&2
    exit 1
}

mkdir -p "$OUT"

# setup.sh must come from this repo (sources CVMFS DELPHI + ROOT 6.34).
# shellcheck disable=SC1091
source "$REPO_ROOT/setup.sh" > /dev/null 2>&1

# -------------------------------------------------------------------- #
# 1. Write a one-line PDL input and run the reader.                     #
# -------------------------------------------------------------------- #
PDL="$OUT/pdlinput_${JOB_ID}"
ROOT_OUT="$OUT/raw_${JOB_ID}.root"
echo "FILE=$AL" > "$PDL"

echo "=== stage 1 : delphi-raw-nanoaod $(basename "$AL") ==="
cd "$OUT"
"$RAW_BIN" --pdlinput "$PDL" --output "$ROOT_OUT" --max-events "$EVENTS" \
    > "$OUT/reader_${JOB_ID}.log" 2>&1
[[ -f "$ROOT_OUT" ]] || { echo "reader failed to produce $ROOT_OUT" >&2; tail -40 "$OUT/reader_${JOB_ID}.log"; exit 1; }
stat -c "  produced: %n  %s bytes" "$ROOT_OUT"
echo "  reader log: $OUT/reader_${JOB_ID}.log"
echo

# -------------------------------------------------------------------- #
# 2. Sanity report via a small ROOT macro.                              #
# -------------------------------------------------------------------- #
MACRO="$OUT/smoke_report.C"
cat > "$MACRO" <<'CPP'
void smoke_report(const char *path) {
    auto r = ROOT::Experimental::RNTupleReader::Open("Events", path);
    Long64_t N = r->GetNEntries();
    std::cout << "entries = " << N << std::endl;
    const char *fields[] = {
        "nTracRaw", "nVdAssocHit", "nVdUnassocHit", "nMtpcRaw",
        "nTrackElement", "nVtx", "nEmShower", "nEmLayer",
        "nHadShower", "nStic", "nMuidRaw", "nElidRaw"
    };
    for (auto f : fields) {
        auto v = r->GetView<std::int16_t>(f);
        long s = 0, mx = 0;
        for (Long64_t i = 0; i < N; ++i) { long x = v(i); s += x; if (x > mx) mx = x; }
        std::cout << "  " << f << ": sum=" << s << " max=" << mx << std::endl;
    }
}
CPP

echo "=== stage 2 : sanity report ==="
root -l -b -q "$MACRO(\"$ROOT_OUT\")" 2>&1 \
    | grep -vE '^Warning|^Processing|startup - Warning'
echo
echo "=== DONE ==="
echo "outputs under $OUT/"
ls -la "$OUT/" | grep -E "${JOB_ID}"
