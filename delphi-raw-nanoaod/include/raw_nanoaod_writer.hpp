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
};

#endif // RAW_NANOAOD_WRITER_HPP
