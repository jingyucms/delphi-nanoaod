# delphi-raw-nanoaod

A PHDST-level reader that converts a DELPHI shortDST (`.sdst`) **or** full-DST
(`.fadana`) into a flat ROOT RNTuple. It is a **sibling** of the SKELANA-based
`delphi-nanoaod` tool in the same repository: they share the `delphi-analysis`
static library, but `delphi-raw-nanoaod` subclasses `phdst::Analysis` instead
of `skelana::Analysis` — so none of SKELANA's aggregation runs, and we walk the
ZEBRA banks directly via `LPHPA` + the `LQ`/`IQ`/`Q` macros.

The point of bypassing SKELANA is to get at **cluster-level and hit-level
information** that SKELANA collapses. For example: SKELANA reports the total
HPC shower energy and a 10-bit layer-hit pattern per track; this reader
surfaces the per-layer energies (EmLayer_*), the conversion-vertex locations
(PhotonConv_* lives on the sibling `feature/photon-collection` branch), the
VD silicon hits with their back-pointers to tracks (VdAssocHit_*), the full
set of per-track detector elements on full-DST input (TrackElement_*), and
the B-field, primary vertex covariance, beamspot dimensions, and the raw MTPC
dE/dx fields that a particle-flow reconstructor needs as inputs.

This README describes what's in the RNTuple, how to build, and how to run.
Design notes and the milestone plan live in `docs/PHDST_RAW_NANOAOD_PLAN.md`.

## Quickstart

```sh
# Inside the cmssw/el9 singularity image with /cvmfs mounted
source /cvmfs/delphi.cern.ch/setup.sh
source /cvmfs/sft.cern.ch/lcg/views/LCG_107/x86_64-el9-gcc13-opt/setup.sh
cd build
cmake ..
make -j4 delphi-raw-nanoaod

./delphi-raw-nanoaod/delphi-raw-nanoaod \
    --pdlinput simana.sdst    --output raw_sdst.root
./delphi-raw-nanoaod/delphi-raw-nanoaod \
    --pdlinput simana.fadana  --output raw_fadana.root
```

Both inputs produce the same RNTuple schema; which collections are populated
depends on the input — see the **shortDST ↔ fullDST** section below.

Quick-look ROOT macro (14 validation histograms):

```sh
root -l -q -b 'scripts/make_plots.C("raw_sdst.root", "plots_dir")'
```

## RNTuple schema

All fields live on a single RNTuple named `Events`. Collections are flat
(one row per `(event, object)`) with back-pointers between them.

### Event-level scalars

| field | source | meaning |
|---|---|---|
| `Event_experimentNumber`, `Event_runNumber`, `Event_fileSequenceNumber`, `Event_eventNumber`, `Event_date`, `Event_time`, `Event_fillNumber` | `/PHCIII/` | event header |
| `Event_bFieldTesla`, `Event_bFieldGevCm` | `BPILOT` | solenoid B (≈ 1.231 T) and the derived GeV/cm factor. `1/R [1/cm] = Event_bFieldGevCm / pT [GeV]` |
| `Event_beamSpot{X,Y,Z}`, `Event_beamSpotSigma{X,Y,Z}`, `Event_beamSpotErrorFlag` | `LQ(LDTOP-25)` | beamspot centre (cm) and Gaussian widths (cm). errorFlag = 0 if the bank is present, -1 if missing |

### Charged tracks (`TracRaw_*`)

One row per charged PA-track (`NINT(Q(LMAIN+8)) ≠ 0`), perigee + fit-quality.
Back-pointers from all other per-track collections use the index in this table.

| field | source | meaning |
|---|---|---|
| `nTracRaw` | count | |
| `TracRaw_paIdx` | PV-PA index | 0-based index of the source PA bank |
| `TracRaw_impactRPhi`, `impactZ`, `theta`, `phi`, `invR` | `QTRAC(4..8)` | perigee parameters |
| `TracRaw_weightMatrix[15]` | `QTRAC(9..23)` | 5×5 symmetric weight matrix |
| `TracRaw_chi2NoVD`, `chi2VD`, `ndfNoVD`, `ndfVD`, `chi2VDHits` | `QMAIN` / `IQ(LMAIN+...)` | fit χ² / ndof split by VD-in / VD-out |
| `TracRaw_trackLength`, `firstPointR`, `firstPointZ` | `QMAIN` | track length (cm) and first measured point |
| `TracRaw_detectorsUsed` | `IQ(LPA+2)` | PSCTRA-documented detector-hit bitmask (see docstring caveat) |
| `TracRaw_charge` | `NINT(Q(LMAIN+8))` | +1 / 0 / -1 (DELPHI charge code decoded) |

### Track elements (`TrackElement_*`) — **fullDST only**

Per-detector 3-D track elements: one row per `PA.TETP` / `TEID` / `TEOD` /
`TEFA` / `TEFB` that the track intersects. **Only populated on `.fadana`
input** — shortDST strips these banks.

| field | meaning |
|---|---|
| `TrackElement_tracRawIdx` | index into `TracRaw_*` |
| `TrackElement_subDetector` | 0 = TPC, 1 = ID, 2 = OD, 3 = FCA, 4 = FCB |
| `TrackElement_dataDescriptor` | bitmask: coord-system + measurement code |
| `TrackElement_coord1..3` | X/Y/Z or R/RPhi/Z at reference point (cm) |
| `TrackElement_theta`, `phi` | direction at the reference point |
| `TrackElement_invPOrPt` | 1/P or 1/Pt (cm⁻¹) |

### VD hits (`VdAssocHit_*`, `VdUnassocHit_*`) — **shortDST only**

Per-hit silicon-strip measurements from the `MVDH` bank at `LQ(LDTOP-21)`.
**Only populated on `.sdst` input** — fullDST does not carry `MVDH`.

Associated hits link back to a `TracRaw` row via `VdAssocHit_tracRawIdx`;
unassociated hits are a flat per-event list.

Fields per hit: `module` (int, sign of Z), `localX` / `R` / `RPhi` (cm),
`signalToNoise` (0 if the shortDST uses the 4-word-per-hit packing).

### TPC per-track summary (`MtpcRaw_*`)

From `PA.MTPC(7)`. Present on both fadana and sdst.

| field | meaning |
|---|---|
| `MtpcRaw_dEdx80Max`, `dEdx80Sigma` | 80%-truncated max-amp dE/dx + sigma |
| `MtpcRaw_dEdx65Max`, `dEdx65Sigma` | 65%-truncated variant |
| `MtpcRaw_dEdx80Integrated` | 80%-trunc integrated amplitude |
| `MtpcRaw_packedPadsSectors` | `nPads + 100*sectorIn + 10000*sectorOut` |
| `MtpcRaw_packedWiresHits` | `nWires + 1000*nCloseClusterHits + 1000000*method` |
| `MtpcRaw_zFitChi2` | χ² of the Z-direction fit |

### Calorimetry

#### EM showers + per-layer: `EmShower_*`, `EmLayer_*`
From `PA.EMCA` / `PA.EMNC`. `EmShower_detector` ∈ {9 HPC, 26 FEMC}. `EmLayer_*`
unpacks the 10-bit HPC layer hit pattern into per-layer raw energies;
**sum(layer energies) is NOT guaranteed to equal EmShower_energy** — the
latter is a calibrated cluster energy with shower-containment corrections,
the former are raw per-layer deposits. They agree to <1 % on high-E clean
EM showers and can differ by factors of a few for soft laterally-spread ones.

#### HAD showers + per-layer-hit: `HadShower_*`, `HadHit_*`
From `PA.HCAL` / `PA.HCNC`. HadHit carries `layer`, `nTowers`, `energy`
per active hit (HAC has 4 layers in DELPHI).

#### STIC: `Stic_*`
Forward small-angle EM. One row per track with a `PA.SSTC` extra-module.

### Photon conversions, tier-2 photon ID, UTER, V0
Exposed on the sibling branch `feature/photon-collection` (the one that also
adds the SKELANA-side `Photon_*` collection). Not in this branch directly
since the two feature branches are intended to stand alone.

### Per-track lepton-ID raw: `MuidRaw_*`, `ElidRaw_*`
`PA.MUID` / `PA.ELID`. Tag, loose-refit χ², hit pattern for muons; tag,
γ-conversion tag, refit 3-momentum for electrons.

### Vertices (`Vtx_*`)

From the PV bank chain at `LQ(LDTOP-1)`. Status bits: 1 = dummy, 2 = secondary,
3 = sec-hadronic, 4 = sim.

| field | meaning |
|---|---|
| `Vtx_position` | (X, Y, Z) in cm |
| `Vtx_chi2`, `Vtx_ndf` | fit quality |
| `Vtx_nOutgoing` | vertex multiplicity |
| `Vtx_massCode` | origin-particle mass code |
| `Vtx_statusBits` | see above |
| `Vtx_errXX / XY / YY / XZ / YZ / ZZ` | symmetric 3×3 error matrix |

## shortDST ↔ fullDST

| info source → | `simana.sdst` (shortDST) | `simana.fadana` (fullDST) |
|---|:---:|:---:|
| `Event_*`, `Vtx_*`, `TracRaw_*` | ✓ | ✓ |
| `EmShower_*` + `EmLayer_*` | ✓ | ✓ |
| `HadShower_*` + `HadHit_*` | ✓ | ✓ |
| `Stic_*`, `MuidRaw_*`, `ElidRaw_*`, `MtpcRaw_*` | ✓ | ✓ |
| `VdAssocHit_*`, `VdUnassocHit_*` (MVDH bank) | **✓** | — |
| `TrackElement_*` (per-detector 3-D hits) | — | **✓** |

For a proper hit-level refit run the reader on **both** files and merge the
two ROOT outputs by `(Event_runNumber, Event_eventNumber)`.

`Delphi-Sim-Pipeline` branch `feature/cmssw-el9-base` preserves both files
in the output directory so the merge is always possible.

## Build notes

Against the CVMFS software stack (LCG_107 + ROOT 6.34 + DELPHI release) in
a `cmssw/el9` singularity image, the target compiles and links cleanly.
See `.github/workflows/build.yml` for the exact recipe — that's what PR CI
uses.

The executable links to libphdst / libdstana / libpxdst / libunfied etc.
plus the CERNLIB set from DELPHI CVMFS, and at runtime the DELPHI binaries
need libgfortran-5, libXm, libXp, libquadmath from the host — the
Delphi-Sim-Pipeline's `container/run_singularity.sh` bind-mounts host
`/lib64` at `/host_lib64` for exactly this reason.

## Scripts

| script | what |
|---|---|
| `scripts/inspect_tracking.C` | aggregate tracking statistics: χ²/ndf, VD hits per track, VD \|R\| histogram, TPC dE/dx, straight-line RPhi residual demo |
| `scripts/make_plots.C` | 14 physics-validation PNGs covering all major collections |

Both are vanilla ROOT C++ macros — run with `root -l -q -b '<name>("input.root", ...)'`.
