// mlpf_export: convert (raw delphi-raw-nanoaod + improved-reco PV refit)
// into a per-event ML-PF style RNTuple, mirroring the CLIC ML-PF dataset
// schema from Zenodo:15062717 (Pata 2024).
//
// Per event we emit three flat tensors:
//
//   X        (n_elements, 17)   per-input-element features (tracks + clusters)
//   ytarget  (n_targets, 13)    truth particles to predict
//   ycand    (n_cands, 13)      legacy reco candidates (for baseline comparison)
//
// CLIC uses the Hungarian-style training where the transformer predicts a
// particle list given X; the loss is computed against ytarget. ycand is the
// legacy PF reconstruction the tagger should beat.
//
// Element-feature schema (matched length-17 to CLIC for compatibility):
//
//   X[i, 0]   element type    1=charged track, 2=EM cluster, 3=HAD cluster, 4=STIC
//   X[i, 1]   pT (track) or E (cluster)        GeV
//   X[i, 2]   eta = -ln(tan(theta/2))
//   X[i, 3]   phi
//   X[i, 4]   charge   (+/-1 for tracks, 0 for clusters)
//   X[i, 5]   d0       (cm, tracks only; 0 for clusters)
//   X[i, 6]   z0       (cm, tracks only; 0 for clusters)
//   X[i, 7]   sigma_d0 (cm, tracks only)
//   X[i, 8]   sigma_z0 (cm, tracks only)
//   X[i, 9]   chi2NormBest (tracks only)
//   X[i, 10]  detector code   (9 = HPC, 26 = FEMC, 0 = other)
//   X[i, 11]  cluster_n_subElements   (e.g. n layers for EM)
//   X[i, 12]  shower energy spread / n_hits
//   X[i, 13]  IP-significance vs PV (tracks; 0 for clusters)
//   X[i, 14]  dz vs PV
//   X[i, 15]  PV.x  (event-level passthrough, repeated per element)
//   X[i, 16]  PV.z  (event-level passthrough)
//
// Particle schema (length-13 for both ytarget and ycand):
//
//   y[i, 0]   PDG id
//   y[i, 1]   charge
//   y[i, 2]   pT
//   y[i, 3]   eta
//   y[i, 4]   phi
//   y[i, 5]   E
//   y[i, 6]   mass
//   y[i, 7]   px
//   y[i, 8]   py
//   y[i, 9]   pz
//   y[i, 10]  vertex.x
//   y[i, 11]  vertex.y
//   y[i, 12]  vertex.z
//
// Trivial v1 of "legacy PF" (ycand):
//   - one charged-pion candidate per preselected charged track,
//   - one photon candidate per EM cluster,
//   - one neutral-hadron candidate per HAD cluster.
// (No track-cluster matching, no particle-flow merging.)
//
// Truth (ytarget) is the GenPart_* status==1 list filtered to charged
// hadrons / leptons / photons / neutral hadrons (status code interpretation
// per JETSET; 1 = stable final-state).

#include <cmath>
#include <cstdint>
#include <iostream>
#include <memory>
#include <vector>
#include <array>

#include <ROOT/RNTuple.hxx>
#include <ROOT/RNTupleModel.hxx>
#include <ROOT/RNTupleReader.hxx>
#include <ROOT/RNTupleWriter.hxx>
#include <Math/Vector3D.h>
#include <Math/Vector4D.h>

using XYZVectorF  = ROOT::Math::DisplacementVector3D<
    ROOT::Math::Cartesian3D<float>,
    ROOT::Math::DefaultCoordinateSystemTag>;
using XYZTVectorF = ROOT::Math::LorentzVector<ROOT::Math::PxPyPzE4D<float>>;
using RNTupleReader = ROOT::Experimental::RNTupleReader;
using RNTupleWriter = ROOT::Experimental::RNTupleWriter;
using RNTupleModel  = ROOT::Experimental::RNTupleModel;

template <typename T>
static std::shared_ptr<T> mk(RNTupleModel &m, const std::string &n,
                             const std::string &d)
{ return m.MakeField<T>({n, d}); }

constexpr int kFeatX = 17;
constexpr int kFeatY = 13;

static float etaFromTheta(float theta) {
    if (theta <= 0 || theta >= M_PI) return 0.f;
    return -std::log(std::tan(0.5f * theta));
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        std::cerr << "Usage: " << argv[0]
                  << " <raw_nanoaod.root> <pv_refit.root> <out_mlpf.root>"
                  << std::endl;
        return 1;
    }
    auto raw = RNTupleReader::Open("Events", argv[1]);
    auto pv  = RNTupleReader::Open("Events", argv[2]);
    const Long64_t N = raw->GetNEntries();
    if (pv->GetNEntries() != N) {
        std::cerr << "Mismatched entry count" << std::endl;
        return 2;
    }

    // Raw nanoaod fields.
    auto vRun  = raw->GetView<int>("Event_runNumber");
    auto vEvt  = raw->GetView<int>("Event_eventNumber");
    auto vIsMC = raw->GetView<std::int8_t>("Event_isMC");
    auto vBgcm = raw->GetView<float>("Event_bFieldGevCm");

    auto vNTracRaw = raw->GetView<std::int16_t>("nTracRaw");
    auto vTracD0   = raw->GetView<std::vector<float>>("TracRaw_impactRPhi");
    auto vTracZ0   = raw->GetView<std::vector<float>>("TracRaw_impactZ");
    auto vTracTh   = raw->GetView<std::vector<float>>("TracRaw_theta");
    auto vTracPhi  = raw->GetView<std::vector<float>>("TracRaw_phi");
    auto vTracInvR = raw->GetView<std::vector<float>>("TracRaw_invR");
    auto vTracQ    = raw->GetView<std::vector<std::int8_t>>("TracRaw_charge");
    auto vTracChi  = raw->GetView<std::vector<float>>("TracRaw_chi2VD");
    auto vTracNdf  = raw->GetView<std::vector<std::int16_t>>("TracRaw_ndfVD");

    auto vNEm   = raw->GetView<std::int16_t>("nEmShower");
    auto vEmDet = raw->GetView<std::vector<std::int8_t>>("EmShower_detector");
    auto vEmE   = raw->GetView<std::vector<float>>("EmShower_energy");
    auto vEmPos = raw->GetView<std::vector<XYZVectorF>>("EmShower_position");
    auto vEmNL  = raw->GetView<std::vector<std::int16_t>>("EmShower_nLayers");

    auto vNHad  = raw->GetView<std::int16_t>("nHadShower");
    auto vHadE  = raw->GetView<std::vector<float>>("HadShower_energy");
    auto vHadDir= raw->GetView<std::vector<XYZVectorF>>("HadShower_direction");
    auto vHadNH = raw->GetView<std::vector<std::int16_t>>("HadShower_nHits");

    auto vNStic = raw->GetView<std::int16_t>("nStic");
    auto vSticE = raw->GetView<std::vector<float>>("Stic_energyFromMain");
    auto vSticD = raw->GetView<std::vector<XYZVectorF>>("Stic_directionFromMain");

    // GenPart_* (truth).
    auto vNGen     = raw->GetView<std::int32_t>("nGenPart");
    auto vGenStat  = raw->GetView<std::vector<std::int16_t>>("GenPart_status");
    auto vGenPdg   = raw->GetView<std::vector<std::int32_t>>("GenPart_pdgId");
    auto vGenP4    = raw->GetView<std::vector<XYZTVectorF>>("GenPart_fourMomentum");
    auto vGenMass  = raw->GetView<std::vector<float>>("GenPart_mass");
    auto vGenVtx   = raw->GetView<std::vector<XYZVectorF>>("GenPart_vertex");

    // PV refit.
    auto vPVR = pv->GetView<XYZVectorF>("PVRefit_position");

    auto model = RNTupleModel::Create();
    auto e_run  = mk<int>(*model,         "Event_runNumber",   "passthrough");
    auto e_evt  = mk<int>(*model,         "Event_eventNumber", "passthrough");
    auto e_isMC = mk<std::int8_t>(*model, "Event_isMC",        "passthrough");
    auto e_pv   = mk<XYZVectorF>(*model,  "PV",                "Refit primary vertex (cm)");
    auto e_X    = mk<std::vector<std::array<float, kFeatX>>>(*model,
            "X",       "Per-element input features (tracks + clusters), 17-D each");
    auto e_yt   = mk<std::vector<std::array<float, kFeatY>>>(*model,
            "ytarget", "Truth particles to predict, 13-D each");
    auto e_yc   = mk<std::vector<std::array<float, kFeatY>>>(*model,
            "ycand",   "Legacy PF candidates (baseline), 13-D each");

    auto writer = RNTupleWriter::Recreate(std::move(model), "Events", argv[3]);

    long sumX = 0, sumYt = 0, sumYc = 0;
    for (Long64_t i = 0; i < N; ++i) {
        *e_run  = vRun(i);
        *e_evt  = vEvt(i);
        *e_isMC = vIsMC(i);
        const auto &PV = vPVR(i);
        *e_pv  = PV;

        e_X->clear();
        e_yt->clear();
        e_yc->clear();

        // ----- X: tracks -----
        const float Bgcm = vBgcm(i);
        const auto &tD0 = vTracD0(i), &tZ0 = vTracZ0(i), &tTh = vTracTh(i);
        const auto &tPh = vTracPhi(i), &tIR = vTracInvR(i);
        const auto &tQ  = vTracQ(i);
        const auto &tCh = vTracChi(i);
        const auto &tNd = vTracNdf(i);
        for (int t = 0; t < vNTracRaw(i); ++t) {
            if (tQ[t] == 0) continue;     // skip neutrals (clusters cover those)
            std::array<float, kFeatX> row{};
            const float pt = (std::fabs(tIR[t]) > 0) ? Bgcm / std::fabs(tIR[t]) : 0.f;
            const float eta = etaFromTheta(tTh[t]);
            const float chiNorm = (tNd[t] > 0) ? tCh[t] / tNd[t] : -1.f;
            row[0]  = 1.f;                  // type = charged track
            row[1]  = pt;
            row[2]  = eta;
            row[3]  = tPh[t];
            row[4]  = static_cast<float>(tQ[t]);
            row[5]  = tD0[t];
            row[6]  = tZ0[t];
            row[7]  = 0.f;   // sigma_d0 not available in raw nanoaod here
            row[8]  = 0.f;
            row[9]  = chiNorm;
            row[10] = 0.f;
            row[11] = 0.f;
            row[12] = 0.f;
            row[13] = 0.f;   // SIP: would need PV here; we skip in trivial v1
            row[14] = tZ0[t] - PV.Z();
            row[15] = PV.X();
            row[16] = PV.Z();
            e_X->push_back(row);

            // Legacy PF: charged pion at this track's kinematics.
            std::array<float, kFeatY> yc{};
            const float p   = pt / std::max(std::sin(tTh[t]), 1e-3f);
            const float pz  = p * std::cos(tTh[t]);
            const float px  = pt * std::cos(tPh[t]);
            const float py  = pt * std::sin(tPh[t]);
            const float Epi = std::sqrt(p*p + 0.13957f*0.13957f);
            yc[0]  = 211.f * tQ[t];     // pi+/- pdg
            yc[1]  = static_cast<float>(tQ[t]);
            yc[2]  = pt;
            yc[3]  = eta;
            yc[4]  = tPh[t];
            yc[5]  = Epi;
            yc[6]  = 0.13957f;
            yc[7]  = px;
            yc[8]  = py;
            yc[9]  = pz;
            // Vertex unknown for legacy PF; use PV.
            yc[10] = PV.X();
            yc[11] = PV.Y();
            yc[12] = PV.Z();
            e_yc->push_back(yc);
        }

        // ----- X: EM clusters -----
        const auto &eDet = vEmDet(i);
        const auto &eE   = vEmE(i);
        const auto &ePos = vEmPos(i);
        const auto &eNL  = vEmNL(i);
        for (int c = 0; c < vNEm(i); ++c) {
            std::array<float, kFeatX> row{};
            const float r = std::sqrt(ePos[c].X()*ePos[c].X() +
                                      ePos[c].Y()*ePos[c].Y() +
                                      ePos[c].Z()*ePos[c].Z());
            const float theta = (r > 0) ? std::acos(ePos[c].Z()/r) : 0.f;
            const float phi   = std::atan2(ePos[c].Y(), ePos[c].X());
            row[0] = 2.f;
            row[1] = eE[c];
            row[2] = etaFromTheta(theta);
            row[3] = phi;
            row[4] = 0.f;
            row[5] = 0.f; row[6] = 0.f; row[7] = 0.f; row[8] = 0.f; row[9] = 0.f;
            row[10] = static_cast<float>(eDet[c]);
            row[11] = static_cast<float>(eNL[c]);
            row[12] = 0.f;
            row[13] = 0.f; row[14] = 0.f;
            row[15] = PV.X(); row[16] = PV.Z();
            e_X->push_back(row);

            std::array<float, kFeatY> yc{};
            yc[0] = 22.f;           // photon
            yc[1] = 0.f;
            yc[2] = eE[c] * std::sin(theta);
            yc[3] = etaFromTheta(theta);
            yc[4] = phi;
            yc[5] = eE[c];
            yc[6] = 0.f;
            yc[7] = eE[c] * std::sin(theta) * std::cos(phi);
            yc[8] = eE[c] * std::sin(theta) * std::sin(phi);
            yc[9] = eE[c] * std::cos(theta);
            yc[10] = PV.X(); yc[11] = PV.Y(); yc[12] = PV.Z();
            e_yc->push_back(yc);
        }

        // ----- X: HAD clusters -----
        const auto &hE   = vHadE(i);
        const auto &hDir = vHadDir(i);
        const auto &hNH  = vHadNH(i);
        for (int c = 0; c < vNHad(i); ++c) {
            std::array<float, kFeatX> row{};
            const float r = std::sqrt(hDir[c].X()*hDir[c].X() +
                                      hDir[c].Y()*hDir[c].Y() +
                                      hDir[c].Z()*hDir[c].Z());
            const float theta = (r > 0) ? std::acos(hDir[c].Z()/r) : 0.f;
            const float phi   = std::atan2(hDir[c].Y(), hDir[c].X());
            row[0] = 3.f;
            row[1] = hE[c];
            row[2] = etaFromTheta(theta);
            row[3] = phi;
            row[10] = 0.f;          // HCAL
            row[11] = static_cast<float>(hNH[c]);
            row[15] = PV.X(); row[16] = PV.Z();
            e_X->push_back(row);

            std::array<float, kFeatY> yc{};
            yc[0] = 130.f;          // K_L (generic neutral hadron)
            yc[1] = 0.f;
            yc[2] = hE[c] * std::sin(theta);
            yc[3] = etaFromTheta(theta);
            yc[4] = phi;
            yc[5] = hE[c];
            yc[6] = 0.4977f;        // K_L mass approx
            yc[7] = hE[c] * std::sin(theta) * std::cos(phi);
            yc[8] = hE[c] * std::sin(theta) * std::sin(phi);
            yc[9] = hE[c] * std::cos(theta);
            yc[10] = PV.X(); yc[11] = PV.Y(); yc[12] = PV.Z();
            e_yc->push_back(yc);
        }

        // ----- X: STIC clusters -----
        const auto &sE   = vSticE(i);
        const auto &sDir = vSticD(i);
        for (int c = 0; c < vNStic(i); ++c) {
            std::array<float, kFeatX> row{};
            const float r = std::sqrt(sDir[c].X()*sDir[c].X() +
                                      sDir[c].Y()*sDir[c].Y() +
                                      sDir[c].Z()*sDir[c].Z());
            const float theta = (r > 0) ? std::acos(sDir[c].Z()/r) : 0.f;
            const float phi   = std::atan2(sDir[c].Y(), sDir[c].X());
            row[0] = 4.f;
            row[1] = sE[c];
            row[2] = etaFromTheta(theta);
            row[3] = phi;
            row[15] = PV.X(); row[16] = PV.Z();
            e_X->push_back(row);
        }

        // ----- ytarget: GenPart status==1 (final-state) -----
        const auto &gSt   = vGenStat(i);
        const auto &gPdg  = vGenPdg(i);
        const auto &gMass = vGenMass(i);
        const auto &gP4   = vGenP4(i);
        const auto &gVtx  = vGenVtx(i);
        for (size_t k = 0; k < gSt.size(); ++k) {
            if (gSt[k] != 1) continue;
            // skip neutrinos (invisible)
            const int apdg = std::abs(gPdg[k]);
            if (apdg == 12 || apdg == 14 || apdg == 16) continue;
            std::array<float, kFeatY> yt{};
            const float p4_pt = std::sqrt(gP4[k].Px()*gP4[k].Px() + gP4[k].Py()*gP4[k].Py());
            const float p4_p  = std::sqrt(p4_pt*p4_pt + gP4[k].Pz()*gP4[k].Pz());
            const float theta = (p4_p > 0) ? std::acos(gP4[k].Pz()/p4_p) : 0.f;
            yt[0] = static_cast<float>(gPdg[k]);
            // charge derived from pdg (very rough; a proper PDG table would be better)
            int q = 0;
            if (apdg == 11 || apdg == 13 || apdg == 15) q = (gPdg[k] > 0) ? -1 : 1;
            else if (apdg == 211 || apdg == 321 || apdg == 2212)
                q = (gPdg[k] > 0) ? 1 : -1;
            yt[1] = q;
            yt[2] = p4_pt;
            yt[3] = etaFromTheta(theta);
            yt[4] = std::atan2(gP4[k].Py(), gP4[k].Px());
            yt[5] = gP4[k].E();
            yt[6] = gMass[k];
            yt[7] = gP4[k].Px();
            yt[8] = gP4[k].Py();
            yt[9] = gP4[k].Pz();
            yt[10] = gVtx[k].X() * 0.1f;   // mm -> cm
            yt[11] = gVtx[k].Y() * 0.1f;
            yt[12] = gVtx[k].Z() * 0.1f;
            e_yt->push_back(yt);
        }

        sumX += e_X->size(); sumYt += e_yt->size(); sumYc += e_yc->size();
        writer->Fill();
    }
    writer.reset();

    std::cout << "Wrote " << argv[3]
              << " : " << N << " events"
              << "  mean X     = " << double(sumX)  / std::max<Long64_t>(N, 1) << " elements/ev"
              << "  mean ytgt  = " << double(sumYt) / std::max<Long64_t>(N, 1) << " particles/ev"
              << "  mean ycand = " << double(sumYc) / std::max<Long64_t>(N, 1) << " candidates/ev"
              << std::endl;
    return 0;
}
