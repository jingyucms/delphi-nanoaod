# delphi-mlpf

ML-PF dataset exporter for DELPHI. Consumes a `delphi-raw-nanoaod` RNTuple
(plus optionally a refit-PV ntuple) and emits a per-event flat-tensor RNTuple
matching the CLIC ML-PF dataset schema (`X` 17-D, `ytarget` 13-D, `ycand`
13-D), so a transformer-based PF reco can train on DELPHI data using the
same recipe as Pata 2024 (zenodo:15062717).

Build with the same recipe as `delphi-raw-nanoaod` (sourced LCG view, then
`g++ -std=c++17 $(root-config --cflags --libs) -lROOTNTuple -lGenVector
mlpf_export.cpp -o mlpf_export`). No DELPHI / SKELANA dependencies.

Run as:
    ./mlpf_export raw_nanoaod.root pv_refit.root out_mlpf.root

Sample numbers on a 100-event Z→bb shortDST: 58 elements/event,
41 truth particles/event, 58 baseline candidates/event.
