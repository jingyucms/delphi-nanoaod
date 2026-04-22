#include "raw_nanoaod_writer.hpp"

#include <cmath>
#include <iostream>

// -----------------------------------------------------------------------------
// Field-registration helpers, copied from delphi-nanoaod/src/nanoaod_writer.cpp.
// Keeping the exact same idiom so the two writers look identical to read.
// -----------------------------------------------------------------------------
template <typename T>
static void MakeField(std::unique_ptr<RNTupleModel> &model,
                      const std::string &name, const std::string &description,
                      std::shared_ptr<T> &ptr)
{
    ptr = model->MakeField<T>({name, description});
}

// -----------------------------------------------------------------------------
// Singleton plumbing. phdst::Analysis::instance_ is the static the user00_()
// shim in phdst_analysis.cpp dispatches through, so the override chain works
// transparently.
// -----------------------------------------------------------------------------
RawNanoAODWriter::RawNanoAODWriter() : output_("raw_nanoaod.root") {}
RawNanoAODWriter::~RawNanoAODWriter() = default;

RawNanoAODWriter *RawNanoAODWriter::getInstance()
{
    if (instance_ == nullptr)
    {
        instance_ = new RawNanoAODWriter();
    }
    return static_cast<RawNanoAODWriter *>(instance_);
}

void RawNanoAODWriter::setOutput(const std::filesystem::path &output)
{
    output_ = output;
}

// -----------------------------------------------------------------------------
// PHDST callbacks.
// -----------------------------------------------------------------------------
void RawNanoAODWriter::user00()
{
    super::user00();

    std::unique_ptr<RNTupleModel> model = RNTupleModel::Create();
    defineEvent(model);
    defineEmShower(model);
    defineHadShower(model);
    defineStic(model);
    defineMuidEl(model);

    writer_ = RNTupleWriter::Recreate(std::move(model), "Events", output_.string());
    std::cout << "RawNanoAODWriter: opened " << output_
              << " (Event, EmShower/Layer, HadShower/Hit, Stic, Muid/ElidRaw)"
              << std::endl;
}

int RawNanoAODWriter::user01()
{
    return super::user01();
}

void RawNanoAODWriter::user02()
{
    super::user02();
    fillEvent();
    fillEmShowers();
    fillHadShowers();
    fillStic();
    fillMuidEl();
    if (!writer_)
    {
        std::cerr << "RawNanoAODWriter::user02: writer_ is null!" << std::endl;
        return;
    }
    writer_->Fill();
    static long filled = 0;
    if (++filled <= 5 || filled % 500 == 0)
    {
        std::cout << "RawNanoAODWriter: filled event " << filled
                  << "  run="    << ph::IIIRUN
                  << "  evt="    << ph::IIIEVT
                  << "  nEmS="   << *nEmShower_
                  << "  nHadS="  << *nHadShower_
                  << "  nStic="  << *nStic_
                  << "  nMuid="  << *nMuidRaw_
                  << "  nElid="  << *nElidRaw_
                  << std::endl;
    }
}

void RawNanoAODWriter::user99()
{
    super::user99();
    if (writer_)
    {
        writer_.reset();   // triggers the RNTuple commit
    }
    std::cout << "RawNanoAODWriter: wrote " << output_ << std::endl;
}

// -----------------------------------------------------------------------------
// Event_* definition + per-event fill (M0 surface).
// -----------------------------------------------------------------------------
void RawNanoAODWriter::defineEvent(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "Event_experimentNumber",   "IIIEXP: experiment number",       Event_experimentNumber_);
    MakeField(model, "Event_runNumber",          "IIIRUN: run number",              Event_runNumber_);
    MakeField(model, "Event_fileSequenceNumber", "IIFILE: file sequence number",    Event_fileSequenceNumber_);
    MakeField(model, "Event_eventNumber",        "IIIEVT: event number",            Event_eventNumber_);
    MakeField(model, "Event_date",               "IIIDAT: event date (yymmdd)",     Event_date_);
    MakeField(model, "Event_time",               "IIITIM: event time (hhmmss)",     Event_time_);
    MakeField(model, "Event_fillNumber",         "IIFILL: LEP fill number",         Event_fillNumber_);
}

void RawNanoAODWriter::fillEvent()
{
    *Event_experimentNumber_   = ph::IIIEXP;
    *Event_runNumber_          = ph::IIIRUN;
    *Event_fileSequenceNumber_ = ph::IIFILE;
    *Event_eventNumber_        = ph::IIIEVT;
    *Event_date_               = ph::IIIDAT;
    *Event_time_               = ph::IIITIM;
    *Event_fillNumber_         = ph::IIFILL;
}

// -----------------------------------------------------------------------------
// EM showers / layers — from the EMNC extra-module under every PA (per-track)
// bank. We walk PV -> PA -> EMNC following the SKELANA PSHEMC pattern
// (skelana.car line 3268).
// -----------------------------------------------------------------------------

void RawNanoAODWriter::defineEmShower(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nEmShower",               "Number of EM-cluster showers (HPC detector code 9, FEMC 26)", nEmShower_);
    MakeField(model, "EmShower_paIdx",          "Index (0-based) of the source PA track bank in this event", EmShower_paIdx_);
    MakeField(model, "EmShower_detector",       "IDET: 9 = HPC (barrel), 26 = EMF/FEMC (endcap merged)",      EmShower_detector_);
    MakeField(model, "EmShower_energy",         "Q(LSHOWR+1) shower energy (GeV)",                            EmShower_energy_);
    MakeField(model, "EmShower_position",       "Q(LSHOWR+2..4) shower impact position (cm)",                 EmShower_position_);
    MakeField(model, "EmShower_nLayers",        "Q(LSHOWR+5) number of HPC layers hit",                       EmShower_nLayers_);
    MakeField(model, "EmShower_layerPattern",   "Q(LSHOWR+6) 10-bit HPC layer hit pattern",                   EmShower_layerPattern_);
    MakeField(model, "EmShower_nLayerEnergies", "Number of EmLayer_* rows emitted for this shower",           EmShower_nLayerEnergies_);

    MakeField(model, "nEmLayer",                "Total EM-layer rows across all showers", nEmLayer_);
    MakeField(model, "EmLayer_emShowerIdx",     "Index into EmShower_* (same row order as the parent)",       EmLayer_emShowerIdx_);
    MakeField(model, "EmLayer_layer",           "HPC layer number 1..NHPLAY",                                 EmLayer_layer_);
    MakeField(model, "EmLayer_energy",
        "Per-layer raw HPC deposit (GeV). NOTE: sum over an EmShower's layer rows "
        "is NOT guaranteed to equal EmShower_energy: the latter is the calibrated, "
        "shower-containment-corrected cluster energy, while the former are raw "
        "per-layer readings. For high-energy EM showers the two agree to <<1%, but "
        "they can differ by factors of a few for soft / laterally-spread showers.",
        EmLayer_energy_);
}

void RawNanoAODWriter::fillEmShowers()
{
    EmShower_paIdx_->clear();
    EmShower_detector_->clear();
    EmShower_energy_->clear();
    EmShower_position_->clear();
    EmShower_nLayers_->clear();
    EmShower_layerPattern_->clear();
    EmShower_nLayerEnergies_->clear();

    EmLayer_emShowerIdx_->clear();
    EmLayer_layer_->clear();
    EmLayer_energy_->clear();

    // NHPLAY = number of HPC layers (9 in DELPHI, verified via pschpc.hpp
    // LENHPC=8 having layer hit pattern at index 6 etc.). Use 10 as a
    // generous loop bound to cover historical 9- and 10-layer variants.
    const int kMaxLayersToScan = 10;

    // Walk: LDTOP -> PV bank(s) -> PA bank(s) -> EMNC.
    if (ph::LDTOP <= 0) { *nEmShower_ = 0; *nEmLayer_ = 0; return; }
    int lpv = ph::LQ(ph::LDTOP - 1);
    int paIdx = 0;

    for (; lpv > 0; lpv = ph::LQ(lpv))
    {
        // Each PV links to a chain of PA (per-track) banks via LQ(LPV - 1).
        for (int lpa = ph::LQ(lpv - 1); lpa > 0; lpa = ph::LQ(lpa), ++paIdx)
        {
            int lemnc = ph::LPHPA("EMNC", lpa, 0);
            if (lemnc == 0) continue;

            int nsidet = static_cast<int>(std::lround(ph::Q(lemnc + 2)));
            int nshowr = nsidet % 100;
            int idet   = nsidet / 100;
            int lshowr = lemnc + 2;

            for (int ns = 0; ns < nshowr; ++ns)
            {
                float e      = ph::Q(lshowr + 1);
                float px     = ph::Q(lshowr + 2);
                float py     = ph::Q(lshowr + 3);
                float pz     = ph::Q(lshowr + 4);
                int   nlay   = static_cast<int>(std::lround(ph::Q(lshowr + 5)));
                int   patt   = static_cast<int>(std::lround(ph::Q(lshowr + 6)));

                EmShower_paIdx_->push_back(static_cast<std::int16_t>(paIdx));
                EmShower_detector_->push_back(static_cast<std::int8_t>(idet));
                EmShower_energy_->push_back(e);
                EmShower_position_->push_back(XYZVectorF(px, py, pz));
                EmShower_nLayers_->push_back(static_cast<std::int16_t>(nlay));
                EmShower_layerPattern_->push_back(patt);

                const std::int16_t showerIdx =
                    static_cast<std::int16_t>(EmShower_paIdx_->size() - 1);

                // Per-layer energies follow the shower header when IDET == 9
                // (HPC, older PXDST flat packing). For EMF or newer packed
                // HPC we skip the walk here; those packings will be handled
                // in follow-up milestones.
                int nEmittedLayers = 0;
                if (idet == 9 && nlay > 0 && nlay <= kMaxLayersToScan)
                {
                    int filledBits = 0;
                    for (int nl = 1; nl <= kMaxLayersToScan && filledBits < nlay; ++nl)
                    {
                        int bit = (patt >> (nl - 1)) & 1;
                        if (!bit) continue;
                        ++filledBits;
                        // Layer energy is Q(LSHOWR + 6 + filledBits).
                        float eLayer = ph::Q(lshowr + 6 + filledBits);
                        EmLayer_emShowerIdx_->push_back(showerIdx);
                        EmLayer_layer_->push_back(static_cast<std::int8_t>(nl));
                        EmLayer_energy_->push_back(eLayer);
                        ++nEmittedLayers;
                    }
                }
                EmShower_nLayerEnergies_->push_back(
                    static_cast<std::int16_t>(nEmittedLayers));

                // Advance LSHOWR to the next shower within this EMNC block
                // (mirrors the SKELANA PSHEMC pattern: 6 header words + NLAY
                // per-layer energies for HPC; 4 words only for EMF).
                if (idet == 9)
                {
                    lshowr += 6 + nlay;
                }
                else
                {
                    lshowr += 4;
                }
            }
        }
    }

    *nEmShower_ = static_cast<std::int16_t>(EmShower_paIdx_->size());
    *nEmLayer_  = static_cast<std::int16_t>(EmLayer_emShowerIdx_->size());
}

// -----------------------------------------------------------------------------
// Hadronic calorimeter — PA.HCNC (shortDST) per PSHHAC (skelana.car L3442).
// -----------------------------------------------------------------------------
void RawNanoAODWriter::defineHadShower(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nHadShower",          "Number of HCAL showers (PA.HCNC)",                   nHadShower_);
    MakeField(model, "HadShower_paIdx",     "PA-track index in the event (0-based)",              HadShower_paIdx_);
    MakeField(model, "HadShower_energy",    "Q(LSHOWR+1) shower energy (GeV)",                    HadShower_energy_);
    MakeField(model, "HadShower_direction", "Q(LSHOWR+2..4) shower 3-direction",                  HadShower_direction_);
    MakeField(model, "HadShower_nHits",     "Q(LSHOWR+5) number of layer-hit pairs",              HadShower_nHits_);
    MakeField(model, "HadShower_nHitRows",  "Number of HadHit_* rows emitted for this shower",    HadShower_nHitRows_);

    MakeField(model, "nHadHit",             "Total HadHit rows", nHadHit_);
    MakeField(model, "HadHit_hadShowerIdx", "Index into HadShower_*",  HadHit_hadShowerIdx_);
    MakeField(model, "HadHit_layer",        "HAC layer number (1..N_HCAL_LAYERS) from packed Q",  HadHit_layer_);
    MakeField(model, "HadHit_nTowers",      "Number of towers in this layer-hit from packed Q",   HadHit_nTowers_);
    MakeField(model, "HadHit_energy",       "Energy of this layer-hit (GeV)",                     HadHit_energy_);
}

void RawNanoAODWriter::fillHadShowers()
{
    HadShower_paIdx_->clear();
    HadShower_energy_->clear();
    HadShower_direction_->clear();
    HadShower_nHits_->clear();
    HadShower_nHitRows_->clear();
    HadHit_hadShowerIdx_->clear();
    HadHit_layer_->clear();
    HadHit_nTowers_->clear();
    HadHit_energy_->clear();

    if (ph::LDTOP <= 0) { *nHadShower_ = 0; *nHadHit_ = 0; return; }
    int paIdx = 0;
    for (int lpv = ph::LQ(ph::LDTOP - 1); lpv > 0; lpv = ph::LQ(lpv))
    {
        for (int lpa = ph::LQ(lpv - 1); lpa > 0; lpa = ph::LQ(lpa), ++paIdx)
        {
            int lhcnc = ph::LPHPA("HCNC", lpa, 0);
            if (lhcnc == 0) continue;

            int nshowr = static_cast<int>(std::lround(ph::Q(lhcnc + 2)));
            int lshowr = lhcnc + 2;

            for (int ns = 0; ns < nshowr; ++ns)
            {
                float e      = ph::Q(lshowr + 1);
                float px     = ph::Q(lshowr + 2);
                float py     = ph::Q(lshowr + 3);
                float pz     = ph::Q(lshowr + 4);
                int   nhits  = static_cast<int>(std::lround(ph::Q(lshowr + 5)));

                HadShower_paIdx_->push_back(static_cast<std::int16_t>(paIdx));
                HadShower_energy_->push_back(e);
                HadShower_direction_->push_back(XYZVectorF(px, py, pz));
                HadShower_nHits_->push_back(static_cast<std::int16_t>(nhits));

                const std::int16_t showerIdx =
                    static_cast<std::int16_t>(HadShower_paIdx_->size() - 1);

                int emitted = 0;
                // Sanity-cap: each hit costs 2 words; a shower's nhits
                // shouldn't be negative or absurdly large.
                if (nhits > 0 && nhits < 200)
                {
                    for (int nl = 1; nl <= nhits; ++nl)
                    {
                        float hitE  = ph::Q(lshowr + 5 + 2*nl - 1);
                        int   packed = static_cast<int>(std::lround(ph::Q(lshowr + 5 + 2*nl)));
                        int   layer = packed / 1000;
                        int   ntow  = packed % 1000;

                        HadHit_hadShowerIdx_->push_back(showerIdx);
                        HadHit_layer_->push_back(static_cast<std::int8_t>(layer));
                        HadHit_nTowers_->push_back(static_cast<std::int16_t>(ntow));
                        HadHit_energy_->push_back(hitE);
                        ++emitted;
                    }
                }
                HadShower_nHitRows_->push_back(static_cast<std::int16_t>(emitted));

                lshowr += 5 + 2*std::max(0, nhits);
            }
        }
    }
    *nHadShower_ = static_cast<std::int16_t>(HadShower_paIdx_->size());
    *nHadHit_    = static_cast<std::int16_t>(HadHit_hadShowerIdx_->size());
}

// -----------------------------------------------------------------------------
// STIC — one row per track with a PA.SSTC extra-module, plus MAIN fields.
// See PSHSTC (skelana.car L4180). We save the commonly-used fields only.
// -----------------------------------------------------------------------------
void RawNanoAODWriter::defineStic(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nStic",                   "Number of tracks with a STIC shower (PA.SSTC)", nStic_);
    MakeField(model, "Stic_paIdx",              "PA-track index",                                Stic_paIdx_);
    MakeField(model, "Stic_energyFromMain",     "Q(LMAIN+6) energy reported by the MAIN block",  Stic_energyFromMain_);
    MakeField(model, "Stic_directionFromMain",  "Q(LMAIN+3..5) direction reported by MAIN",      Stic_directionFromMain_);
    MakeField(model, "Stic_numHitTowers",       "Q(LSSTC+2) encoded tower-hit count (/10)",      Stic_numHitTowers_);
}

void RawNanoAODWriter::fillStic()
{
    Stic_paIdx_->clear();
    Stic_energyFromMain_->clear();
    Stic_directionFromMain_->clear();
    Stic_numHitTowers_->clear();

    if (ph::LDTOP <= 0) { *nStic_ = 0; return; }
    int paIdx = 0;
    for (int lpv = ph::LQ(ph::LDTOP - 1); lpv > 0; lpv = ph::LQ(lpv))
    {
        for (int lpa = ph::LQ(lpv - 1); lpa > 0; lpa = ph::LQ(lpa), ++paIdx)
        {
            int lsstc = ph::LPHPA("SSTC", lpa, 0);
            if (lsstc == 0) continue;

            int lmain = ph::LPHPA("MAIN", lpa, 0);
            if (lmain == 0) continue;

            // SSTC payload: header + per-shower block (SKELANA only looks at
            // the first shower, with the loop commented out in PSHSTC).
            int lshowr = lsstc + 2;

            Stic_paIdx_->push_back(static_cast<std::int16_t>(paIdx));
            Stic_energyFromMain_->push_back(ph::Q(lmain + 6));
            Stic_directionFromMain_->push_back(
                XYZVectorF(ph::Q(lmain + 3), ph::Q(lmain + 4), ph::Q(lmain + 5)));
            int rawHitCount = static_cast<int>(std::lround(ph::Q(lshowr + 2) / 10.0));
            Stic_numHitTowers_->push_back(static_cast<std::int16_t>(rawHitCount));
        }
    }
    *nStic_ = static_cast<std::int16_t>(Stic_paIdx_->size());
}

// -----------------------------------------------------------------------------
// Per-track Muon + Electron raw ID (PA.MUID + PA.ELID).
// Mirrors PSCMUD / PSCELD but emitted by our PHDST reader directly, without
// going through SKELANA. Useful when the SKELANA-based pipeline is not in
// play but you still want the standard per-track lepton-ID words.
// -----------------------------------------------------------------------------
void RawNanoAODWriter::defineMuidEl(std::unique_ptr<RNTupleModel> &model)
{
    MakeField(model, "nMuidRaw",            "Number of tracks with PA.MUID (muon-ID)", nMuidRaw_);
    MakeField(model, "MuidRaw_paIdx",       "PA-track index",                          MuidRaw_paIdx_);
    MakeField(model, "MuidRaw_tag",         "Q(LMUID+1) NINT: MUCAL2 tag",             MuidRaw_tag_);
    MakeField(model, "MuidRaw_looseChi2",   "Q(LMUID+2) global chi2 of very loose refit", MuidRaw_looseChi2_);
    MakeField(model, "MuidRaw_hitPattern",  "Q(LMUID+3) NINT: hit pattern with inefficiencies", MuidRaw_hitPattern_);

    MakeField(model, "nElidRaw",            "Number of tracks with PA.ELID (electron-ID)", nElidRaw_);
    MakeField(model, "ElidRaw_paIdx",       "PA-track index",                              ElidRaw_paIdx_);
    MakeField(model, "ElidRaw_tag",         "Q(LELID+1) NINT: ELECID tag (0..5)",          ElidRaw_tag_);
    MakeField(model, "ElidRaw_gammaConvTag","Q(LELID+2) NINT: gamma-conversion tag",       ElidRaw_gammaConvTag_);
    MakeField(model, "ElidRaw_refitMomentum","Q(LELID+3..5): electron refit 3-momentum at vertex", ElidRaw_refitMomentum_);
}

void RawNanoAODWriter::fillMuidEl()
{
    MuidRaw_paIdx_->clear();
    MuidRaw_tag_->clear();
    MuidRaw_looseChi2_->clear();
    MuidRaw_hitPattern_->clear();

    ElidRaw_paIdx_->clear();
    ElidRaw_tag_->clear();
    ElidRaw_gammaConvTag_->clear();
    ElidRaw_refitMomentum_->clear();

    if (ph::LDTOP <= 0) { *nMuidRaw_ = 0; *nElidRaw_ = 0; return; }
    int paIdx = 0;
    for (int lpv = ph::LQ(ph::LDTOP - 1); lpv > 0; lpv = ph::LQ(lpv))
    {
        for (int lpa = ph::LQ(lpv - 1); lpa > 0; lpa = ph::LQ(lpa), ++paIdx)
        {
            int lmuid = ph::LPHPA("MUID", lpa, 0);
            if (lmuid > 0)
            {
                MuidRaw_paIdx_->push_back(static_cast<std::int16_t>(paIdx));
                MuidRaw_tag_->push_back(
                    static_cast<std::int32_t>(std::lround(ph::Q(lmuid + 1))));
                MuidRaw_looseChi2_->push_back(ph::Q(lmuid + 2));
                MuidRaw_hitPattern_->push_back(
                    static_cast<std::int32_t>(std::lround(ph::Q(lmuid + 3))));
            }

            int lelid = ph::LPHPA("ELID", lpa, 0);
            if (lelid > 0)
            {
                ElidRaw_paIdx_->push_back(static_cast<std::int16_t>(paIdx));
                ElidRaw_tag_->push_back(
                    static_cast<std::int32_t>(std::lround(ph::Q(lelid + 1))));
                ElidRaw_gammaConvTag_->push_back(
                    static_cast<std::int32_t>(std::lround(ph::Q(lelid + 2))));
                ElidRaw_refitMomentum_->push_back(XYZVectorF(
                    ph::Q(lelid + 3), ph::Q(lelid + 4), ph::Q(lelid + 5)));
            }
        }
    }
    *nMuidRaw_ = static_cast<std::int16_t>(MuidRaw_paIdx_->size());
    *nElidRaw_ = static_cast<std::int16_t>(ElidRaw_paIdx_->size());
}
