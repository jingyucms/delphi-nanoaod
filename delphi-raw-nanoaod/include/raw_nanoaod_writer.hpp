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
};

#endif // RAW_NANOAOD_WRITER_HPP
