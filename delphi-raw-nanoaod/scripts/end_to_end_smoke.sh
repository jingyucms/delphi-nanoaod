#!/bin/bash
# End-to-end smoke test: Pythia8 -> DELSIM -> (shortDST, fullDST) ->
# delphi-raw-nanoaod -> ROOT, plus a brief sanity report.
#
# Intended to run locally inside the cmssw/el9 singularity image with
# /cvmfs visible. Driven by the run_singularity.sh wrapper from
# Delphi-Sim-Pipeline (feature/cmssw-el9-base), which preserves both
# simana_<job>.sdst and simana_<job>.fadana in the output directory.
#
# Usage
# -----
#   end_to_end_smoke.sh [--events N] [--out DIR] [--config PATH]
#                       [--sim-pipeline DIR] [--raw-bin PATH] [--image PATH]
#
# Defaults
# --------
#   --events       100
#   --out          /tmp/phdst_smoke
#   --config       $SIM_PIPELINE/config_z_mm.txt
#   --sim-pipeline $HOME/Codes/ee_Z_tatagamma/Delphi-Sim-Pipeline  (edit to taste)
#   --raw-bin      $(this repo)/build_raw/delphi-raw-nanoaod/delphi-raw-nanoaod
#   --image        $HOME/Codes/ee_Z_tatagamma/images/cmssw-el9.sif
#
# Any non-default value can also come from an env var of the same
# uppercased name, e.g. `EVENTS=300 ./end_to_end_smoke.sh` .

set -eo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO_ROOT="$(cd "$HERE/../.." && pwd)"

EVENTS=${EVENTS:-100}
OUT=${OUT:-/tmp/phdst_smoke}
SIM_PIPELINE=${SIM_PIPELINE:-$HOME/Codes/ee_Z_tatagamma/Delphi-Sim-Pipeline}
CONFIG=${CONFIG:-$SIM_PIPELINE/config_z_mm.txt}
RAW_BIN=${RAW_BIN:-$REPO_ROOT/build_raw/delphi-raw-nanoaod/delphi-raw-nanoaod}
IMAGE=${IMAGE:-$HOME/Codes/ee_Z_tatagamma/images/cmssw-el9.sif}
JOB_ID=${JOB_ID:-smoke_$(date +%Y%m%d_%H%M%S)}

while [[ $# -gt 0 ]]; do
    case "$1" in
        --events)       EVENTS=$2;       shift 2;;
        --out)          OUT=$2;          shift 2;;
        --config)       CONFIG=$2;       shift 2;;
        --sim-pipeline) SIM_PIPELINE=$2; shift 2;;
        --raw-bin)      RAW_BIN=$2;      shift 2;;
        --image)        IMAGE=$2;        shift 2;;
        --job-id)       JOB_ID=$2;       shift 2;;
        -h|--help)      sed -n '2,30p' "$0"; exit 0;;
        *) echo "Unknown flag: $1" >&2; exit 2;;
    esac
done

echo "=== config ==="
echo "  events       = $EVENTS"
echo "  out          = $OUT"
echo "  sim-pipeline = $SIM_PIPELINE"
echo "  config       = $CONFIG"
echo "  raw-bin      = $RAW_BIN"
echo "  image        = $IMAGE"
echo "  job-id       = $JOB_ID"
echo

mkdir -p "$OUT"

# -------------------------------------------------------------------- #
# 1. Pythia8 + DELSIM via the cvmfs-backed singularity wrapper.         #
#    Produces $OUT/simana_$JOB_ID.sdst and .fadana.                     #
# -------------------------------------------------------------------- #
wrap="$SIM_PIPELINE/container/run_singularity.sh"
[[ -x "$wrap" ]] || { echo "missing wrapper: $wrap" >&2; exit 1; }
[[ -x "$RAW_BIN" ]] || { echo "missing raw-reader binary: $RAW_BIN" >&2; exit 1; }
[[ -f "$IMAGE" ]] || { echo "missing singularity image: $IMAGE" >&2; exit 1; }
[[ -f "$CONFIG" ]] || { echo "missing pythia config: $CONFIG" >&2; exit 1; }

export IMAGE_DIR="$(dirname "$IMAGE")"
echo "=== stage 1 : Pythia8 + DELSIM ($EVENTS events) ==="
"$wrap" "$EVENTS" "$JOB_ID" "$OUT" "$CONFIG" | tail -12
echo

SDST="$OUT/simana_${JOB_ID}.sdst"
FADANA="$OUT/simana_${JOB_ID}.fadana"
for f in "$SDST" "$FADANA"; do
    [[ -f "$f" ]] || { echo "missing expected output: $f" >&2; exit 1; }
done
echo "got:"
ls -la "$SDST" "$FADANA"
echo

# -------------------------------------------------------------------- #
# 2. delphi-raw-nanoaod on both inputs.                                 #
#    Needs /cvmfs + host /lib64 for the DELPHI fortran libs at runtime. #
# -------------------------------------------------------------------- #
SDST_ROOT="$OUT/raw_${JOB_ID}_sdst.root"
FADANA_ROOT="$OUT/raw_${JOB_ID}_fadana.root"

run_reader () {
    local input=$1
    local output=$2
    local tag=$3
    echo "=== stage 2 ($tag): delphi-raw-nanoaod $(basename "$input") ==="
    singularity exec \
        --cleanenv \
        --bind /cvmfs \
        --bind /lib64:/host_lib64:ro \
        --bind "$OUT:$OUT" \
        --bind "$(dirname "$RAW_BIN")":/binroot:ro \
        "$IMAGE" bash -c "
            set -eo pipefail
            source /cvmfs/delphi.cern.ch/setup.sh > /dev/null 2>&1
            source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh
            export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/host_lib64
            cd $OUT
            ln -sf $(basename "$input") T.FSEQ1
            /binroot/$(basename "$RAW_BIN") \
                --pdlinput '$input' --output '$output' \
                --max-events $EVENTS > /tmp/raw_${tag}.log 2>&1
        "
    [[ -f "$output" ]] || { echo "reader failed to produce $output" >&2; exit 1; }
    stat -c "  produced: %n  %s bytes" "$output"
}
run_reader "$SDST"   "$SDST_ROOT"   "sdst"
run_reader "$FADANA" "$FADANA_ROOT" "fadana"
echo

# -------------------------------------------------------------------- #
# 3. Sanity-report via a small ROOT macro.                              #
# -------------------------------------------------------------------- #
MACRO=/tmp/smoke_report.C
cat > "$MACRO" <<'CPP'
void smoke_report(const char *path, const char *tag) {
    auto r = ROOT::Experimental::RNTupleReader::Open("Events", path);
    Long64_t N = r->GetNEntries();
    std::cout << "--- " << tag << " (" << path << ") -- entries = " << N << std::endl;
    for (const char *f : {"nTracRaw", "nVdAssocHit", "nMtpcRaw",
                          "nTrackElement", "nVtx", "nEmShower", "nEmLayer",
                          "nHadShower", "nStic", "nMuidRaw", "nElidRaw"}) {
        auto v = r->GetView<std::int16_t>(f);
        long s = 0; for (Long64_t i = 0; i < N; ++i) s += v(i);
        std::cout << "  total " << f << " = " << s << std::endl;
    }
}
CPP

echo "=== stage 3 : sanity report ==="
singularity exec \
    --cleanenv \
    --bind /cvmfs \
    --bind /lib64:/host_lib64:ro \
    --bind "$OUT:$OUT" \
    --bind "$MACRO:$MACRO:ro" \
    "$IMAGE" bash -c "
        set +u
        source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh
        export LD_LIBRARY_PATH=\$LD_LIBRARY_PATH:/host_lib64
        root -l -q -b '
            .L $MACRO
            smoke_report(\"$SDST_ROOT\",   \"shortDST\");
            smoke_report(\"$FADANA_ROOT\", \"fullDST \");
        ' 2>&1 | grep -vE 'Setting LC_|startup - Warning|Processing'
    "
echo
echo "=== DONE ==="
echo "outputs under $OUT/"
ls -la "$OUT/" | grep -E "${JOB_ID}"
