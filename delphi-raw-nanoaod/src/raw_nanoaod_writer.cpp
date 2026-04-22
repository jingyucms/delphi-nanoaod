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

    writer_ = RNTupleWriter::Recreate(std::move(model), "Events", output_.string());
    std::cout << "RawNanoAODWriter: opened " << output_ << " with "
              << "Event_* + EmShower_* + EmLayer_*" << std::endl;
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
                  << "  run="   << ph::IIIRUN
                  << "  evt="   << ph::IIIEVT
                  << "  nEmS="  << *nEmShower_
                  << "  nEmL="  << *nEmLayer_
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
