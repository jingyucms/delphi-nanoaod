// Worked example: what the raw nanoAOD gives you for track-refit work.
//
// Uses only fields from feature/phdst-raw-reader:
//   Event_bFieldTesla, Event_bFieldGevCm
//   TracRaw_{impactRPhi, impactZ, theta, phi, invR, chi2VD, ndfVD, charge}
//   VdAssocHit_{tracRawIdx, R, RPhi, signalToNoise}
//   MtpcRaw_{dEdx80Max, dEdx80Sigma}
//
// For each event, we dump:
//   - B field (should read 1.231 T for DELPHI)
//   - number of charged tracks, VD hits, MTPC rows
//   - for charged tracks, number of VD hits actually assigned to it
//   - aggregated chi2/ndf distribution (peaks around 1 for a healthy fit)
//   - dE/dx distribution (peaks near 1 MIP for minimum-ionising tracks)
//   - VD hit S/N distribution (VD's 5-layer silicon is designed for S/N > 5)
// No refit, no helix propagation. A helix-based residual demo would
// need DELPHI's sign-of-kappa and helix-centering convention pinned
// down against a reference track, which is a separate exercise.
//
// Run:
//   root -l -q -b 'vd_residuals.C("/path/to/raw_nanoaod.root", 200)'

#include "ROOT/RNTuple.hxx"
#include "ROOT/RNTupleReader.hxx"

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <iostream>
#include <vector>

void inspect_tracking(const char *path, Long64_t max_events = 200)
{
    using namespace ROOT::Experimental;
    auto reader = RNTupleReader::Open("Events", path);
    Long64_t N = reader->GetNEntries();
    Long64_t lim = std::min<Long64_t>(max_events, N);
    std::cout << "entries   : " << N << std::endl;
    std::cout << "processing: " << lim << " events" << std::endl;

    auto vB    = reader->GetView<float>("Event_bFieldTesla");
    auto vBgc  = reader->GetView<float>("Event_bFieldGevCm");
    auto vNTR  = reader->GetView<std::int16_t>("nTracRaw");
    auto vNVD  = reader->GetView<std::int16_t>("nVdAssocHit");
    auto vNMT  = reader->GetView<std::int16_t>("nMtpcRaw");

    auto vC2   = reader->GetView<std::vector<float>>("TracRaw_chi2VD");
    auto vNd   = reader->GetView<std::vector<std::int16_t>>("TracRaw_ndfVD");
    auto vTI   = reader->GetView<std::vector<std::int16_t>>("VdAssocHit_tracRawIdx");
    auto vSN   = reader->GetView<std::vector<float>>("VdAssocHit_signalToNoise");
    auto vR    = reader->GetView<std::vector<float>>("VdAssocHit_R");
    auto vDedx = reader->GetView<std::vector<float>>("MtpcRaw_dEdx80Max");

    std::cout << "B field    : " << vB(0) << " T,  BGeVCm = " << vBgc(0) << std::endl;

    long tot_tr = 0, tot_hit = 0, tot_mt = 0;
    long tot_good_chi2 = 0, tot_fit = 0;
    double sum_chi2_over_ndf = 0.0;

    // per-track VD hit count histogram
    std::array<long, 15> nHitHist{}; // 0..13 plus overflow
    // S/N histogram (buckets at 0,2,4,8,16,32, +inf)
    const double snEdges[] = {0.0, 2.0, 4.0, 8.0, 16.0, 32.0};
    std::array<long, 7> snHist{};
    // dE/dx histogram in MIPs (buckets 0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0, +inf)
    const double eEdges[] = {0.0, 0.5, 1.0, 1.5, 2.0, 3.0, 5.0};
    std::array<long, 8> eHist{};
    // VD hit R histogram (DELPHI layer radii 6.3, 9.0, 10.9 cm)
    const double rEdges[] = {0.0, 5.0, 7.0, 9.5, 12.0, 20.0};
    std::array<long, 7> rHist{};

    for (Long64_t i = 0; i < lim; ++i)
    {
        tot_tr  += vNTR(i);
        tot_hit += vNVD(i);
        tot_mt  += vNMT(i);

        const auto &chi = vC2(i);
        const auto &nd  = vNd(i);
        for (std::size_t k = 0; k < chi.size(); ++k)
        {
            if (nd[k] > 0 && chi[k] > 0)
            {
                double x = chi[k] / double(nd[k]);
                sum_chi2_over_ndf += x;
                if (x < 5.0) ++tot_good_chi2;
                ++tot_fit;
            }
        }

        // VD hits per track — build a reverse index.
        const auto &ti = vTI(i);
        std::vector<int> hitsPerTrack(vNTR(i), 0);
        for (auto t : ti) if (t >= 0 && t < (int)hitsPerTrack.size()) hitsPerTrack[t]++;
        for (int h : hitsPerTrack)
        {
            int b = std::min<int>(h, (int)nHitHist.size() - 1);
            ++nHitHist[b];
        }

        const auto &sn = vSN(i);
        for (auto s : sn)
        {
            int b = 6;
            for (int e = 0; e < 6; ++e) if (s < snEdges[e+1 < 6 ? e+1 : 5] && e < 5) { b = e; break; }
            // Simpler classification:
            if      (s < snEdges[1]) b = 0;
            else if (s < snEdges[2]) b = 1;
            else if (s < snEdges[3]) b = 2;
            else if (s < snEdges[4]) b = 3;
            else if (s < snEdges[5]) b = 4;
            else                      b = 5;
            ++snHist[b];
        }

        const auto &r = vR(i);
        for (auto rr : r)
        {
            double ar = std::fabs(rr);
            int b = 6;
            if      (ar < rEdges[1]) b = 0;
            else if (ar < rEdges[2]) b = 1;
            else if (ar < rEdges[3]) b = 2;
            else if (ar < rEdges[4]) b = 3;
            else if (ar < rEdges[5]) b = 4;
            else                     b = 5;
            ++rHist[b];
        }

        const auto &e = vDedx(i);
        for (auto ev : e)
        {
            int b = 7;
            if (ev > 0)
            {
                if      (ev < eEdges[1]) b = 0;
                else if (ev < eEdges[2]) b = 1;
                else if (ev < eEdges[3]) b = 2;
                else if (ev < eEdges[4]) b = 3;
                else if (ev < eEdges[5]) b = 4;
                else if (ev < eEdges[6]) b = 5;
                else                     b = 6;
            }
            ++eHist[b];
        }
    }

    std::cout << std::endl;
    std::cout << "average per event:" << std::endl;
    std::cout << "  charged tracks (TracRaw) : " << double(tot_tr)/lim << std::endl;
    std::cout << "  VD associated hits       : " << double(tot_hit)/lim << std::endl;
    std::cout << "  MTPC rows                : " << double(tot_mt)/lim << std::endl;
    std::cout << std::endl;
    std::cout << "charged-track fits with VD:" << std::endl;
    std::cout << "  tracks with ndfVD>0      : " << tot_fit << std::endl;
    std::cout << "  of which chi2/ndf < 5    : " << tot_good_chi2
              << " (" << 100.0*tot_good_chi2/std::max(1L,tot_fit) << " %)" << std::endl;
    if (tot_fit) std::cout
        << "  <chi2/ndf>               : " << sum_chi2_over_ndf/tot_fit << std::endl;

    std::cout << std::endl << "VD hits per track:" << std::endl;
    for (std::size_t k = 0; k < nHitHist.size(); ++k)
        if (nHitHist[k] > 0) std::cout << "  " << k
            << (k+1 == nHitHist.size() ? "+" : " ")
            << " hits: " << nHitHist[k] << " tracks" << std::endl;

    std::cout << std::endl << "VD hit |R| (cm, DELPHI layers at ~6.3 / 9.0 / 10.9):" << std::endl;
    const char *rLabels[] = {"<5", "5-7", "7-9.5", "9.5-12", "12-20", ">=20"};
    for (int k = 0; k < 6; ++k)
        std::cout << "  " << rLabels[k] << ": " << rHist[k] << std::endl;

    std::cout << std::endl << "VD hit signal-to-noise:" << std::endl;
    const char *snLabels[] = {"<2", "2-4", "4-8", "8-16", "16-32", ">=32"};
    for (int k = 0; k < 6; ++k)
        std::cout << "  " << snLabels[k] << ": " << snHist[k] << std::endl;

    std::cout << std::endl << "MtpcRaw dE/dx (80% trunc, MIP-normalised):" << std::endl;
    const char *eLabels[] = {"0-0.5", "0.5-1", "1-1.5", "1.5-2", "2-3", "3-5", ">5", "==0 (no dE/dx)"};
    for (int k = 0; k < 8; ++k)
        std::cout << "  " << eLabels[k] << ": " << eHist[k] << std::endl;

    // -----------------------------------------------------------------
    // Straight-line RPhi residual: predicted RPhi at cylinder r_hit is
    //   RPhi_pred = r_hit * phi0 + d0
    // to first order in d0/r_hit and ignoring curvature. For a 1 GeV track
    // in the DELPHI 1.23 T field, the neglected curvature term is
    // O(r_hit^2 / 2 R_curve) ≈ 2 mm at r_hit = 11 cm. So good residuals
    // look like a ~ mm-scale Gaussian dominated by curvature plus a ~ μm
    // tail of the true VD resolution. A full-helix prediction would close
    // that last factor, but needs the DELPHI sign-of-kappa convention
    // pinned down; left as a follow-up.
    // -----------------------------------------------------------------
    auto vD0_tr  = reader->GetView<std::vector<float>>("TracRaw_impactRPhi");
    auto vPhi_tr = reader->GetView<std::vector<float>>("TracRaw_phi");
    auto vTI2    = reader->GetView<std::vector<std::int16_t>>("VdAssocHit_tracRawIdx");
    auto vHR2    = reader->GetView<std::vector<float>>("VdAssocHit_R");
    auto vHP2    = reader->GetView<std::vector<float>>("VdAssocHit_RPhi");

    std::vector<double> residuals;
    for (Long64_t i = 0; i < lim; ++i) {
        const auto &d0s = vD0_tr(i);
        const auto &phs = vPhi_tr(i);
        const auto &ti  = vTI2(i);
        const auto &hR  = vHR2(i);
        const auto &hP  = vHP2(i);
        for (std::size_t k = 0; k < ti.size(); ++k) {
            int t = ti[k];
            if (t < 0 || t >= (int)d0s.size()) continue;
            double r = hR[k];
            if (r <= 0) continue;          // Z-measuring hits have R < 0
            double pred = r * phs[t] + d0s[t];
            double res  = pred - hP[k];
            double wrap = 2.0 * M_PI * r;
            while (res >  0.5 * wrap) res -= wrap;
            while (res < -0.5 * wrap) res += wrap;
            residuals.push_back(res);
        }
    }
    if (!residuals.empty()) {
        double mean = 0.0;
        for (double r : residuals) mean += r;
        mean /= residuals.size();
        double var = 0.0;
        for (double r : residuals) var += (r - mean) * (r - mean);
        var /= residuals.size();
        double rms = std::sqrt(var);
        long w5mm = 0, w2mm = 0, w02mm = 0;
        for (double r : residuals) {
            if (std::abs(r) < 0.5)  ++w5mm;
            if (std::abs(r) < 0.2)  ++w2mm;
            if (std::abs(r) < 0.02) ++w02mm;
        }
        std::cout << std::endl << "Straight-line RPhi residuals (pred = r*phi0 + d0):" << std::endl;
        std::cout << "  hits used    : " << residuals.size() << std::endl;
        std::cout << "  mean         : " << mean*1e4 << " um" << std::endl;
        std::cout << "  RMS          : " << rms *1e4 << " um" << std::endl;
        std::cout << "  within 5 mm  : " << w5mm  << " ("
                  << 100.0*w5mm /residuals.size() << " %)" << std::endl;
        std::cout << "  within 2 mm  : " << w2mm  << " ("
                  << 100.0*w2mm /residuals.size() << " %)" << std::endl;
        std::cout << "  within 0.2 mm: " << w02mm << " ("
                  << 100.0*w02mm/residuals.size() << " %)" << std::endl;
        std::cout << "  note: O(r^2/2R) curvature term ~ 2 mm at 1 GeV, "
                  << "7 mm at 300 MeV -- a full helix would tighten RMS." << std::endl;
    }
}
