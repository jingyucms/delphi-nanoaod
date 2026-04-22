#ifndef RAW_NANOAOD_WRITER_HPP
#define RAW_NANOAOD_WRITER_HPP

// PHDST-level raw-hit nanoAOD writer.
// Subclasses phdst::Analysis (NOT skelana::Analysis) so that no SKELANA
// aggregation runs: we read the ZEBRA banks directly. M0 just surfaces the
// PHCIII event-header fields; M1+ add PA.EMCA / PA.HCAL / PA.STIC / … hit
// collections. See docs/PHDST_RAW_NANOAOD_PLAN.md.

#include <filesystem>
#include <memory>
#include <vector>

#include "phdst_analysis.hpp"
#include "phdst.hpp"

#include <array>

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <Math/Vector3D.h>

using XYZVectorF = ROOT::Math::DisplacementVector3D<
    ROOT::Math::Cartesian3D<float>,
    ROOT::Math::DefaultCoordinateSystemTag>;

using RNTupleWriter = ROOT::Experimental::RNTupleWriter;
using RNTupleModel = ROOT::Experimental::RNTupleModel;

namespace ph = phdst;

class RawNanoAODWriter : public ph::Analysis
{
private:
    typedef ph::Analysis super;

public:
    RawNanoAODWriter(const RawNanoAODWriter &) = delete;
    RawNanoAODWriter &operator=(const RawNanoAODWriter &) = delete;
    virtual ~RawNanoAODWriter();
    static RawNanoAODWriter *getInstance();
    void setOutput(const std::filesystem::path &output);

protected:
    RawNanoAODWriter();
    virtual void user00() override;
    virtual int  user01() override;
    virtual void user02() override;
    virtual void user99() override;

private:
    void defineEvent(std::unique_ptr<RNTupleModel> &model);
    void fillEvent();
    void defineEmShower(std::unique_ptr<RNTupleModel> &model);
    void fillEmShowers();   // walks PA banks, emits EmShower_* + EmLayer_*
    void defineHadShower(std::unique_ptr<RNTupleModel> &model);
    void fillHadShowers();  // PA.HCNC -> HadShower_* + HadHit_*
    void defineStic(std::unique_ptr<RNTupleModel> &model);
    void fillStic();        // PA.SSTC -> Stic_* (one row per track with STIC)
    void defineMuidEl(std::unique_ptr<RNTupleModel> &model);
    void fillMuidEl();      // PA.MUID + PA.ELID -> MuidRaw_* + ElidRaw_*
    void defineTrac(std::unique_ptr<RNTupleModel> &model);
    void fillTrac();        // PA.TRAC + PA.MAIN -> TracRaw_* per charged track

    std::filesystem::path              output_;
    std::unique_ptr<RNTupleWriter>     writer_;

    // Event_* fields (M0): from PHCIII common block.
    std::shared_ptr<int>   Event_experimentNumber_;
    std::shared_ptr<int>   Event_runNumber_;
    std::shared_ptr<int>   Event_fileSequenceNumber_;
    std::shared_ptr<int>   Event_eventNumber_;
    std::shared_ptr<int>   Event_date_;
    std::shared_ptr<int>   Event_time_;
    std::shared_ptr<int>   Event_fillNumber_;

    // EM-shower collections (M2).
    //
    // Sourced from the EMNC extra-module under every PA (per-track) bank in
    // the shortDST. Each EMNC may contain 0, 1, ... showers; each shower
    // unpacks an energy, a position-vector, a 10-bit HPC layer pattern, and
    // optionally (for HPC showers in older PXDST versions) one per-layer
    // energy per set bit.
    std::shared_ptr<std::int16_t>                         nEmShower_;
    std::shared_ptr<std::vector<std::int16_t>>            EmShower_paIdx_;         // index of source PA track (0-based)
    std::shared_ptr<std::vector<std::int8_t>>             EmShower_detector_;      // 9 = HPC, 26 = EMF(FEMC)
    std::shared_ptr<std::vector<float>>                   EmShower_energy_;        // Q(LSHOWR+1)
    std::shared_ptr<std::vector<XYZVectorF>>              EmShower_position_;      // Q(LSHOWR+2..4) in cm
    std::shared_ptr<std::vector<std::int16_t>>            EmShower_nLayers_;       // Q(LSHOWR+5)
    std::shared_ptr<std::vector<std::int32_t>>            EmShower_layerPattern_;  // Q(LSHOWR+6) 10-bit HPC hit pattern
    std::shared_ptr<std::vector<std::int16_t>>            EmShower_nLayerEnergies_;// number of EmLayer rows for this shower

    std::shared_ptr<std::int16_t>                         nEmLayer_;
    std::shared_ptr<std::vector<std::int16_t>>            EmLayer_emShowerIdx_;    // index into EmShower_*
    std::shared_ptr<std::vector<std::int8_t>>             EmLayer_layer_;          // 1..NHPLAY (DELPHI convention)
    std::shared_ptr<std::vector<float>>                   EmLayer_energy_;         // per-layer deposit

    // --- Hadronic-calorimeter collections (M3): PA.HCNC shortDST extra-module.
    // Layout per PSHHAC (skelana.car line 3442): per shower,
    //   Q(LSHOWR+1) = shower E
    //   Q(LSHOWR+2..4) = 3-momentum (x, y, z) of the shower direction
    //   Q(LSHOWR+5) = NLHITS = number of layer-hit pairs
    //   Q(LSHOWR+5+2*nl-1) = energy of hit nl
    //   Q(LSHOWR+5+2*nl)   = packed 1000*layer + nTowers
    // Advance by 5 + 2*NLHITS.
    std::shared_ptr<std::int16_t>                         nHadShower_;
    std::shared_ptr<std::vector<std::int16_t>>            HadShower_paIdx_;
    std::shared_ptr<std::vector<float>>                   HadShower_energy_;
    std::shared_ptr<std::vector<XYZVectorF>>              HadShower_direction_;
    std::shared_ptr<std::vector<std::int16_t>>            HadShower_nHits_;
    std::shared_ptr<std::vector<std::int16_t>>            HadShower_nHitRows_;

    std::shared_ptr<std::int16_t>                         nHadHit_;
    std::shared_ptr<std::vector<std::int16_t>>            HadHit_hadShowerIdx_;
    std::shared_ptr<std::vector<std::int8_t>>             HadHit_layer_;     // 1..N_HCAL_LAYERS
    std::shared_ptr<std::vector<std::int16_t>>            HadHit_nTowers_;
    std::shared_ptr<std::vector<float>>                   HadHit_energy_;

    // --- STIC small-angle EM collection (M4): PA.SSTC shortDST.
    // One row per track with a STIC shower (SKELANA's CCC DO NS loop is
    // commented out so there's effectively one shower per track). No layer
    // decomposition — STIC is thin enough that SKELANA only keeps totals.
    std::shared_ptr<std::int16_t>                         nStic_;
    std::shared_ptr<std::vector<std::int16_t>>            Stic_paIdx_;
    std::shared_ptr<std::vector<float>>                   Stic_energyFromMain_;   // Q(LMAIN+6)
    std::shared_ptr<std::vector<XYZVectorF>>              Stic_directionFromMain_;
    std::shared_ptr<std::vector<std::int16_t>>            Stic_numHitTowers_;

    // --- Muon / Electron per-track raw ID (M5): PA.MUID + PA.ELID records.
    // Keep it simple — one row per track that has that extra-module present,
    // carrying the raw ID words. The Muid_* / Elid_* collections in the
    // SKELANA-based delphi-nanoaod carry the parsed/refitted versions; these
    // sit alongside them and expose the un-interpreted PSCMUD/PSCELD words
    // for downstream ML / calibration work.
    std::shared_ptr<std::int16_t>                         nMuidRaw_;
    std::shared_ptr<std::vector<std::int16_t>>            MuidRaw_paIdx_;
    std::shared_ptr<std::vector<std::int32_t>>            MuidRaw_tag_;        // Q(LMUID+1) NINT = tag
    std::shared_ptr<std::vector<float>>                   MuidRaw_looseChi2_;  // Q(LMUID+2)
    std::shared_ptr<std::vector<std::int32_t>>            MuidRaw_hitPattern_; // Q(LMUID+3) NINT

    std::shared_ptr<std::int16_t>                         nElidRaw_;
    std::shared_ptr<std::vector<std::int16_t>>            ElidRaw_paIdx_;
    std::shared_ptr<std::vector<std::int32_t>>            ElidRaw_tag_;             // Q(LELID+1) NINT
    std::shared_ptr<std::vector<std::int32_t>>            ElidRaw_gammaConvTag_;    // Q(LELID+2) NINT
    std::shared_ptr<std::vector<XYZVectorF>>              ElidRaw_refitMomentum_;   // Q(LELID+3..5)

    // --- Track raw (M6): PA.TRAC + PA.MAIN. One row per charged track
    // (Q(LMAIN+8) != 0). Follows SKELANA PSHTRA / PSCTRA (stdcdes.car
    // +KEEP,PSCTRA). The 20-word TRAC payload carries the perigee
    // parameters and the 15-element symmetric weight matrix.
    std::shared_ptr<std::int16_t>                         nTracRaw_;
    std::shared_ptr<std::vector<std::int16_t>>            TracRaw_paIdx_;
    std::shared_ptr<std::vector<float>>                   TracRaw_impactRPhi_;    // QTRAC( 4)
    std::shared_ptr<std::vector<float>>                   TracRaw_impactZ_;       // QTRAC( 5)
    std::shared_ptr<std::vector<float>>                   TracRaw_theta_;         // QTRAC( 6)
    std::shared_ptr<std::vector<float>>                   TracRaw_phi_;           // QTRAC( 7)
    std::shared_ptr<std::vector<float>>                   TracRaw_invR_;          // QTRAC( 8)  curvature w/ sign at perigee
    std::shared_ptr<std::vector<std::array<float, 15>>>   TracRaw_weightMatrix_;  // QTRAC(9..23) 5x5 symmetric
    std::shared_ptr<std::vector<float>>                   TracRaw_trackLength_;   // |Q(LMAIN+9)| (sign flipped by SKELANA)
    std::shared_ptr<std::vector<std::int32_t>>            TracRaw_detectorsUsed_; // IQ(LPA+2) bits 1-VD,2-ID,3-TPC,4-OD,5-FCA,6-FCB
    std::shared_ptr<std::vector<float>>                   TracRaw_firstPointR_;   // |Q(LMAIN+23..24)|
    std::shared_ptr<std::vector<float>>                   TracRaw_firstPointZ_;   // Q(LMAIN+25)
    std::shared_ptr<std::vector<float>>                   TracRaw_chi2NoVD_;      // Q(LMAIN+16)
    std::shared_ptr<std::vector<float>>                   TracRaw_chi2VD_;        // Q(LMAIN+26)
    std::shared_ptr<std::vector<std::int16_t>>            TracRaw_ndfNoVD_;       // Q(LMAIN+17)
    std::shared_ptr<std::vector<std::int16_t>>            TracRaw_ndfVD_;         // Q(LMAIN+27)
    std::shared_ptr<std::vector<float>>                   TracRaw_chi2VDHits_;    // Q(LMAIN+18)
    std::shared_ptr<std::vector<std::int8_t>>             TracRaw_charge_;        // sign of Q(LMAIN+8)
};

#endif // RAW_NANOAOD_WRITER_HPP
