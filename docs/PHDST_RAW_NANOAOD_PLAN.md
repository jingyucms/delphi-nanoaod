# Plan: a PHDST-level raw nanoAOD reader (`feature/phdst-raw-reader`)

## Context

`delphi-nanoaod` currently sits **on top of SKELANA** — `NanoAODWriter` inherits
`skelana::Analysis`, SKELANA fills high-level COMMON blocks (VECP, PSCHPC,
QPHOT, QPI0, ...), and we marshal those into the RNTuple. This is ergonomic and
covers the standard DELPHI physics objects (tracks, reconstructed showers, V0s,
π⁰ tags, …) but SKELANA **aggregates** calorimeter data: it exposes a shower's
*total* energy and a 9- or 10-bit *layer hit pattern*, not the per-layer (or
per-cluster / per-cell) energies that particle-flow and ML-based reconstruction
want.

The shortDST already **carries** the finer-grained data — see
`/cvmfs/delphi.cern.ch/attic/cc7-gfortran-no-f2c/simana/v94c/dat/shortdst.des`
block `PA.EMCA(2)`, lines 1155-1247: every EM shower has a per-layer or
per-cluster tail with energies and positions (the precise packing depends on
the PXDST version). The information is accessible from C++ via the **PHDST**
layer that already ships inside the repo (`delphi-analysis/include/phdst/`) —
we just need to bypass SKELANA and walk the ZEBRA banks directly.

This plan is for a second, sibling executable (`delphi-raw-nanoaod`) that
reads the SDST at the PHDST level and writes a sibling RNTuple with raw
calorimeter / tracker hits. It is deliberately *not* a replacement for the
SKELANA-based `delphi-nanoaod`; it is a complement for analyses that need
what SKELANA throws away. The two can run on the same SDST and produce
separate ROOT files which downstream code merges by run/event number.

## What's already reusable in the repo

`delphi-analysis/include/phdst/` has the C++ surface for talking to PHDST:

| file | what it gives us |
|---|---|
| `phdst/functions.hpp` | `PHDST`, `PHSET`, `IPHPIC`, `PHRTY`, `TIMED`, `TIMEX` — steering |
| `phdst/uxcom.hpp` | `LQ(i)`, `IQ(i)`, `Q(i)` — the classic ZEBRA pointer-chase macros into the main store |
| `phdst/uxlink.hpp` | `LRTOP`, `LSTOP`, `LTTOP`, `LBKTOP`, … — the top-of-store link words |
| `phdst/phciii.hpp` | `IIIEXP / IIIRUN / IIIEVT / IIIDAT / IIITIM / IIFILL` — event header |
| `phdst/phgen.hpp`, `phdst/pxchdr.hpp` | miscellaneous |

`delphi-analysis/src/phdst_analysis.cpp` wraps the
`user00/user01/user02/user99` callbacks as a `phdst::Analysis` class, exactly
the same user-override pattern SKELANA uses — except no SKELANA commons are
filled. Subclassing `phdst::Analysis` instead of `skelana::Analysis` is the
clean "bypass SKELANA" path.

## What we do NOT yet have and need to add

1. **Bank-access bindings.** PHDST's Fortran side has routines to look up a
   bank by name and return a ZEBRA link word (e.g. `LMDGET`, `MZLINK`,
   `IPHPIL`, `LPHPA`). None are exposed in `phdst/functions.hpp` today. Add
   thin `extern "C"` wrappers there, modelled on the existing six.
2. **A new `delphi-analysis/include/phdst/banks.hpp`** (or similar) that
   documents and names the extra-module layout constants we actually walk
   into — `PA.EMCA = 2`, `PA.HCAL = 3`, `PA.STIC = 19`, `PA.SSTC = 33`,
   `PA.TRAC = 8`, `PA.OTRK = 31`, `PA.PHOT = 30`, etc., matching the
   `#TITLE` lines in `shortdst.des`.
3. **The new executable**, `delphi-raw-nanoaod`, under
   `delphi-nanoaod/delphi-raw-nanoaod/` — new CMake target sibling to the
   existing `delphi-nanoaod/delphi-nanoaod/`. Rough layout:
   ```
   delphi-nanoaod/delphi-raw-nanoaod/
     CMakeLists.txt
     include/raw_nanoaod_writer.hpp
     src/delphi-raw-nanoaod.cpp
     src/raw_nanoaod_writer.cpp
     src/skelana.cra.in           # still needed so nypatchy pulls in STDCDES
                                  # for the few Fortran-side bits we reuse
   ```
   `RawNanoAODWriter` subclasses `phdst::Analysis`. `user00` defines the
   RNTuple model; `user02` walks the PA extra-module tree for each event
   via the bank-access helpers and pushes rows into the per-collection
   vectors; `user99` commits the RNTuple.

## Proposed RNTuple schema

Per-event collections, each one an "unrolled" view of the matching shortDST
bank. Names follow the `Camel_snakeSuffix` convention used everywhere else.

| collection | source (shortdst.des title) | per-row fields |
|---|---|---|
| `Event_*` | PILOT.IDEN / PHCIII | run, event, date, time, exp, fill, DST version, PXDST version |
| `EmShower_*` | `PA.EMCA(2)` loop over showers | detector (9=HPC, 26=FEMC), E, x, y, z, θ, φ, massId, nLayers, layerPattern |
| `EmLayer_*` | `PA.EMCA(2)` inner DO NLAY_HPC loop (PXDST ≤ 262) | emShowerIdx, layer #, E, nChannels |
| `EmCluster_*` | `PA.EMCA(2)` inner DO NCLU loop (PXDST ≥ 263) | emShowerIdx, packed/unpacked (layer, x, y, z, width_in_z, E) |
| `HadShower_*` | `PA.HCAL(3)` | detector (10=HAC), E, x, y, z, θ, φ, massId |
| `HadLayer_*` | `PA.HCAL(3)` inner layer loop | hadShowerIdx, layer #, E, nTowers |
| `SticShower_*` | `PA.STIC(19)` / `PA.SSTC(33)` | side, tower, energy, x, y, z |
| `Track_*` | `PA.TRAC(8)`, augmented by `PA.OTRK(31)` | partial-track segments with per-detector hits |
| `Mu_*` | `PA.MU(4)` | muon-chamber hit row per track |
| `El_*` | `PA.EL(5)` | electron-ID row per track |
| `Tpc_*` | `PA.MTPC(7)` | TPC segment info |

Branch field types: `std::int16_t` or `std::int32_t` for counters / indices,
`float` for energies and coordinates, `bool` for flags. We reuse the same
`MakeField`/`fillVector` helpers from `nanoaod_writer.cpp`.

**Shape**: one entry per event (like `delphi-nanoaod`), so downstream code
joins the raw ntuple with the SKELANA ntuple by `(Event_runNumber,
Event_eventNumber)`. No attempt to reuse the *same* RNTuple as the SKELANA
path — keeping the schemas separate avoids entangling the two reader classes.

## Concrete milestones

Rough order, each one an independently shippable commit:

**M0 — scaffold.** New CMake target, empty `RawNanoAODWriter` subclassing
`phdst::Analysis`, writes just the PHCIII-based `Event_*` fields. Tests:
runs on a Z→μμ SDST, output has nEntries matching the SKELANA nanoAOD, and
`Event_runNumber` matches.

**M1 — bank-access helpers.** Add `phdst/banks.hpp` + the `LMDGET` /
`IPHPIL` bindings in `phdst/functions.hpp`. Prove it works by logging the
lengths of `PA.EMCA` / `PA.HCAL` / `PA.STIC` per event.

**M2 — EmShower + EmLayer/EmCluster.** Walk `PA.EMCA(2)` per event. For
each shower emit one `EmShower_*` row (detector / energy / position / layer
pattern) and, per inner layer or cluster, one `EmLayer_*` or `EmCluster_*`
row indexed back to the parent via `emShowerIdx`. Handle the three
sub-formats (PXDST < 253, 253 ≤ v < 263, v ≥ 263); start with v ≥ 263
(modern data) and add the older variants only when the reader encounters
them.
Validation: `sum(EmLayer_E) for emShowerIdx == k` ≈ `EmShower_E[k]` to float
precision, and per-event sum of `EmShower_E` on HPC showers matches the
SKELANA-side `Part_hpcTotalShowerEnergy` aggregated over neutrals.

**M3 — HadShower + HadLayer** from `PA.HCAL(3)`. Same pattern.

**M4 — SticShower / SticTower** from `PA.STIC(19)` (barrel SAT-like) and
`PA.SSTC(33)` (small-angle).

**M5 — Track / OTrk / TPC.** Walk `PA.TRAC(8)` for the track-element bank
and `PA.OTRK(31)` for outer-tracker hits; optionally `PA.MTPC(7)`. This is
the heaviest and can be deferred behind a CMake / config flag.

**M6 — Mu + El.** `PA.MU(4)`, `PA.EL(5)`: raw per-layer muon-chamber and
electron-ID information that SKELANA collapses to a single tag.

Everything above milestone M1 is gated on the bank-access helpers landing
and being tested. M2 alone is enough to claim "per-layer HPC output ready"
for particle-flow use.

## fadana vs shortDST — complementary information

DELPHI's reconstruction pipeline generates two artefacts per job:

* **`simana.fadana`** — DELANA's full-DST output. Carries the per-track
  3-D track elements in `PA.TETP` (TPC), `PA.TEID` (ID), `PA.TEOD` (OD),
  `PA.TEFA` (FCA), `PA.TEFB` (FCB). This is what `TrackElement_*` reads.
  It does **not** carry the `MVDH` VD-hit bank at `LQ(LDTOP-21)`.
* **`simana.sdst`** — `shortdst.exe`'s compressed shortDST. Carries
  `MVDH` (consumed into `VdAssocHit_*` / `VdUnassocHit_*`). Strips the
  `PA.TE*` banks entirely.

Both files have the same PHDST on-disk layout, so the same
`delphi-raw-nanoaod` binary reads either one — the populated collections
depend on which input you feed it. For a proper hit-based refit run the
reader twice and merge the two ROOT outputs by
`(Event_runNumber, Event_eventNumber)`.

`Delphi-Sim-Pipeline` branch `feature/cmssw-el9-base` has been updated
to preserve `simana.fadana` alongside `simana.sdst` in every job's
output directory, so this merge is always possible.

## Out of scope (explicitly)

- **True raw RDST banks.** Only the shortDST (processed) banks are in our
  reach — SDST is already calibrated / clustered. Pulling genuine raw
  readout data (unclustered HPC cells etc.) would require a separate
  pipeline reading the raw data tape, which we do not have access to.
- **Re-running reconstruction.** If M2 shows unacceptable cluster quality
  we do not re-do the clustering; we expose what's there.
- **Calorimeter geometry.** We save cluster `(x, y, z)` as reported; we do
  not carry full cell-geometry look-ups. Downstream can map position →
  cell via an external geometry file if needed.
- **Merging with the SKELANA-based RNTuple at write time.** Kept as two
  separate files; join downstream.

## Risks / open questions

1. **PXDST version spread.** A single SDST usually comes from one PXDST
   version, but the dataset as a whole spans v<253 through v≥270. The
   reader must detect the version (stored in `PILOT.VERS`) and branch on
   it. This is all in `shortdst.des` but will take iteration to get right.
2. **DST type.** MiniDST / shortDST / fullDST have overlapping but distinct
   bank trees. We target **shortDST** first (matches the `delphi-nanoaod`
   input), but leave a `dstType` field in `Event_*` so mini or full can be
   added later.
3. **ZEBRA bank-name case / padding.** PHDST expects 4-char bank names; our
   bindings need to pad or truncate correctly, which the existing `PHDST`
   wrapper does by `std::strncpy`. The same pattern should work for
   `LMDGET` and friends, but it's worth a unit test on startup.
4. **Output size.** A Z→hadrons event has O(10^2) HPC clusters × O(10)
   layers or clusters = ~10^3 rows per event per sub-detector. At 16 bytes
   per row and 17 k PHDST-respooled entries per SDST, a 200-event SDST
   could blow up to a few hundred MB of raw-nanoAOD. Feasible but not
   free; include a `max_hit_pattern_bits` / `skip_tracks` config flag.
5. **The `T.FSEQ1` fallback in PHDST** (the same one documented in
   `Delphi-Sim-Pipeline/.github/workflows/smoke-test.yml`) respools each
   input SDST ~30x, producing duplicate events in the output tree. Either
   fix PHDST (hard), detect and skip duplicate (run, event) pairs (the
   shipping workaround), or document that the raw-nanoAOD inherits the
   duplication and downstream must dedup.

## Verification / acceptance tests

Once M2 is done, rerun the end-to-end hadronic smoke (the same 200-event
`Z→light quarks` SDST used for the Photon tier-2 verification at
`/tmp/smoke_had/simana_hadronic.sdst`). Cross-checks:

1. `Event_runNumber` / `Event_eventNumber` join 1:1 between the skelana
   and raw ntuples (modulo the PHDST respooling — same duplication factor
   on both sides).
2. Per neutral `Part_*` of mass-code `2` (photon-ish) with
   `Part_hpcTotalShowerEnergy > 0`, there should be at least one matching
   `EmShower_*` row on HPC with the same `(theta, phi)` and an energy
   agreeing within a few %.
3. `sum_over_layers(EmLayer_E)` equals `EmShower_E` to float precision for
   modern-PXDST events.

Until these three checks pass, the reader is pre-alpha.

## Estimated effort

Rough, honest numbers:

- M0: < 1 day (scaffold + one scalar branch).
- M1: 1-2 days (LMDGET binding, find-bank plumbing).
- M2: 3-5 days (PXDST version handling is the real cost).
- M3 / M4: 1-2 days each, once M2 works.
- M5 / M6: 2-3 days each.
- Overall "useful for particle flow on HPC": ~1 week of focused work to get
  M0-M2; another 1-2 weeks for the rest.

## Status (updated)

Branch `feature/phdst-raw-reader` is off `main` (rebased to be independent
of `feature/photon-collection` per the user's ask), and all planned
milestones have landed. Summary:

| milestone | status | what it added |
|---|---|---|
| M0 | ✓ `e0f98c7` | scaffold, `Event_*` from PHCIII |
| M1 + M2 | ✓ `f66a5af` | `LPHPA` binding, UXLINK alignment fix, `EmShower_* + EmLayer_*` |
| M3 + M4 + M5 | ✓ `9d6f4d9` | `HadShower_* + HadHit_*`, `Stic_*`, `MuidRaw_*`, `ElidRaw_*` |
| M6 | ✓ `b94ce00` | `TracRaw_*` per-track perigee + 5×5 weight matrix |
| build fix | ✓ `a674a70` | `pscvda.hpp` + `pscvdu.hpp` + FASTJETDIR default |
| CI | ✓ `3d3ebba` | `.github/workflows/build.yml` |
| M7 | ✓ `441316c` | VD hit collections, `TrackElement_*` skeleton, `MtpcRaw_*` |
| M8 | ✓ `2834b31` | `Event_bField*` + inspection macro |
| straight-line residual demo | ✓ `1ca5926` | first-order RPhi residual in `inspect_tracking.C` |
| `scripts/make_plots.C` | ✓ `b42d835` | 8 physics-validation PNGs |
| M9 | ✓ `8fbb8bb` | beamspot + `Vtx_*` + 6 more plots |
| fadana / shortDST docs | ✓ `8b94eae` | confirm `TrackElement_*` populates on `.fadana` input |
| `README.md` + this status | in this commit | |

For a hit-level refit, run the reader twice (on `.sdst` and `.fadana`) and
merge by `(Event_runNumber, Event_eventNumber)`. Every other collection that
matters for a DELPHI analysis is exposed.
