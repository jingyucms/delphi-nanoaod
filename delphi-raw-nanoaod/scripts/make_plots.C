// Produce a small set of physics-validation PNGs from a raw nanoAOD file.
// All plots are made with vanilla TH1/TCanvas so they render the same in
// any ROOT 6.x install.
//
// Run:
//   root -l -q -b 'make_plots.C("/path/to/raw_nanoaod.root", "/tmp/plots")'

#include "ROOT/RNTuple.hxx"
#include "ROOT/RNTupleReader.hxx"
#include "Math/Vector3D.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TStyle.h"
#include "TLine.h"
#include "TLatex.h"

using XYZVectorF = ROOT::Math::DisplacementVector3D<
    ROOT::Math::Cartesian3D<float>,
    ROOT::Math::DefaultCoordinateSystemTag>;

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <filesystem>
#include <iostream>
#include <string>
#include <vector>

static void save(TCanvas &c, const std::string &dir, const std::string &stem) {
    std::filesystem::create_directories(dir);
    std::string path = dir + "/" + stem + ".png";
    c.SaveAs(path.c_str());
    std::cout << "wrote " << path << std::endl;
}

void make_plots(const char *path, const char *outdir = "/tmp/raw_plots",
                Long64_t max_events = 200)
{
    using namespace ROOT::Experimental;
    gStyle->SetOptStat(1110);
    gStyle->SetPadLeftMargin(0.14);
    gStyle->SetPadBottomMargin(0.13);

    auto r = RNTupleReader::Open("Events", path);
    Long64_t N = r->GetNEntries();
    Long64_t lim = std::min(max_events, N);
    std::cout << "entries = " << N << "  processing " << lim << std::endl;

    // Views
    auto vBT   = r->GetView<float>("Event_bFieldTesla");
    auto vNTR  = r->GetView<std::int16_t>("nTracRaw");
    auto vNVD  = r->GetView<std::int16_t>("nVdAssocHit");
    auto vC2   = r->GetView<std::vector<float>>("TracRaw_chi2VD");
    auto vNd   = r->GetView<std::vector<std::int16_t>>("TracRaw_ndfVD");
    auto vInvR = r->GetView<std::vector<float>>("TracRaw_invR");
    auto vTh   = r->GetView<std::vector<float>>("TracRaw_theta");
    auto vPhi  = r->GetView<std::vector<float>>("TracRaw_phi");
    auto vD0   = r->GetView<std::vector<float>>("TracRaw_impactRPhi");
    auto vTI   = r->GetView<std::vector<std::int16_t>>("VdAssocHit_tracRawIdx");
    auto vHR   = r->GetView<std::vector<float>>("VdAssocHit_R");
    auto vHP   = r->GetView<std::vector<float>>("VdAssocHit_RPhi");
    auto vDEdx = r->GetView<std::vector<float>>("MtpcRaw_dEdx80Max");
    auto vEmE  = r->GetView<std::vector<float>>("EmShower_energy");
    auto vHadE = r->GetView<std::vector<float>>("HadShower_energy");

    // Histograms
    TH1F hNtrac ("hNtrac",  ";charged tracks / event;events",       40, -0.5, 59.5);
    TH1F hVDR   ("hVDR",    ";VD hit |R| (cm);hits",                120, 4.0, 13.0);
    TH1F hChi2  ("hChi2",   ";#chi^{2}/ndf (with VD);tracks",        80, 0.0, 8.0);
    TH1F hPt    ("hPt",     ";p_{T} [GeV/c];tracks",                 50, 0.0, 20.0);
    TH1F hDEdx  ("hDEdx",   ";TPC dE/dx (80% trunc, MIP);tracks",    80, 0.0, 5.0);
    TH1F hEm    ("hEm",     ";EmShower energy [GeV];showers",        60, 0.0, 40.0);
    TH1F hHad   ("hHad",    ";HadShower energy [GeV];showers",       60, 0.0, 20.0);
    TH1F hRes   ("hRes",    ";straight-line RPhi residual [mm];hits",
                                                               100, -30.0, 30.0);

    auto bfield = vBT(0);   // 1.231 T
    auto bGevCm = r->GetView<float>("Event_bFieldGevCm")(0);

    for (Long64_t i = 0; i < lim; ++i) {
        hNtrac.Fill(vNTR(i));
        const auto &chi = vC2(i); const auto &nd = vNd(i);
        const auto &kap = vInvR(i); const auto &th = vTh(i);
        const auto &phi = vPhi(i); const auto &d0 = vD0(i);
        for (std::size_t k = 0; k < kap.size(); ++k) {
            if (nd[k] > 0 && chi[k] > 0) hChi2.Fill(chi[k] / nd[k]);
            if (std::abs(kap[k]) > 0) {
                // pT [GeV] = BGEVCM / |kappa [1/cm]|
                double pT = bGevCm / std::abs(kap[k]);
                if (pT < 50.0) hPt.Fill(pT);
            }
        }
        const auto &hR = vHR(i);
        for (auto rr : hR) hVDR.Fill(std::abs(rr));
        const auto &ti = vTI(i); const auto &hP = vHP(i);
        for (std::size_t k = 0; k < ti.size(); ++k) {
            int t = ti[k];
            if (t < 0 || t >= (int)kap.size()) continue;
            double rh = hR[k];
            if (rh <= 0) continue;
            double pred = rh * phi[t] + d0[t];
            double res  = pred - hP[k];
            double wrap = 2.0 * M_PI * rh;
            while (res >  0.5*wrap) res -= wrap;
            while (res < -0.5*wrap) res += wrap;
            hRes.Fill(res * 10.0);  // cm -> mm
        }
        for (auto d : vDEdx(i)) if (d > 0) hDEdx.Fill(d);
        for (auto e : vEmE(i))  hEm.Fill(e);
        for (auto e : vHadE(i)) hHad.Fill(e);
    }

    std::string out = outdir;

    { TCanvas c("c","",800,600); hNtrac.Draw(); save(c, out, "01_tracks_per_event"); }
    { TCanvas c("c","",800,600); hVDR.SetFillColor(kAzure-9); hVDR.Draw();
      TLatex txt; txt.SetNDC(); txt.SetTextSize(0.032);
      txt.DrawLatex(0.56, 0.80, "DELPHI VD layers:");
      txt.DrawLatex(0.56, 0.76, "  r = 6.3, 9.0, 10.9 cm");
      save(c, out, "02_vd_hit_R"); }
    { TCanvas c("c","",800,600); hChi2.Draw(); save(c, out, "03_chi2_per_ndf"); }
    { TCanvas c("c","",800,600); c.SetLogy(); hPt.Draw();
      TLatex txt; txt.SetNDC(); txt.SetTextSize(0.032);
      char buf[128];
      snprintf(buf, sizeof(buf), "B = %.3f T (from BPILOT)", bfield);
      txt.DrawLatex(0.55, 0.82, buf);
      save(c, out, "04_pT_spectrum"); }
    { TCanvas c("c","",800,600); hDEdx.SetFillColor(kOrange-2); hDEdx.Draw();
      save(c, out, "05_dEdx"); }
    { TCanvas c("c","",800,600); c.SetLogy(); hEm.Draw();
      save(c, out, "06_em_shower_energy"); }
    { TCanvas c("c","",800,600); c.SetLogy(); hHad.Draw();
      save(c, out, "07_had_shower_energy"); }
    { TCanvas c("c","",800,600); hRes.Draw();
      TLatex txt; txt.SetNDC(); txt.SetTextSize(0.032);
      txt.DrawLatex(0.55, 0.82, "First-order only:");
      txt.DrawLatex(0.55, 0.78, "  RPhi_{pred} = r #phi_{0} + d_{0}");
      txt.DrawLatex(0.55, 0.72, "curvature term");
      txt.DrawLatex(0.55, 0.68, "  O(r^{2}/2R_{c}) ~ 2 mm @ 1 GeV");
      save(c, out, "08_vd_rphi_residual"); }

    // -----------------------------------------------------------------
    // Extra: vertex + beamspot plots (M9 collections).
    // -----------------------------------------------------------------
    auto vNV   = r->GetView<std::int16_t>("nVtx");
    auto vVpos = r->GetView<std::vector<XYZVectorF>>("Vtx_position");
    auto vVnO  = r->GetView<std::vector<std::int16_t>>("Vtx_nOutgoing");
    auto vVst  = r->GetView<std::vector<std::int32_t>>("Vtx_statusBits");
    auto vBX   = r->GetView<float>("Event_beamSpotX");
    auto vBY   = r->GetView<float>("Event_beamSpotY");
    auto vBZ   = r->GetView<float>("Event_beamSpotZ");

    TH1F hNv    ("hNv",    ";vertices / event;events",                   30, -0.5, 29.5);
    TH1F hVz    ("hVz",    ";vertex Z (cm);vertices",                    80, -40.0, 40.0);
    TH1F hVrho  ("hVrho",  ";vertex #rho = #sqrt{x^{2}+y^{2}} (cm);vertices",
                                                                          60, 0.0, 6.0);
    TH1F hVprim ("hVprim", ";vertex #rho (cm, primary only);vertices",    60, 0.0, 0.2);
    TH1F hVno   ("hVno",   ";vertex multiplicity (nOutgoing);vertices",   40, 0.5, 40.5);
    TH1F hBx    ("hBx",    ";beam spot X per event (mm);events",          80, -2.0, 2.0);
    TH1F hBy    ("hBy",    ";beam spot Y per event (mm);events",          80, -0.5, 0.5);
    TH1F hBz    ("hBz",    ";beam spot Z per event (mm);events",          80, -5.0, 5.0);

    for (Long64_t i = 0; i < lim; ++i) {
        hNv.Fill(vNV(i));
        hBx.Fill(vBX(i) * 10.0);
        hBy.Fill(vBY(i) * 10.0);
        hBz.Fill(vBZ(i) * 10.0);
        const auto &pos = vVpos(i);
        const auto &nO  = vVnO(i);
        const auto &st  = vVst(i);
        for (std::size_t k = 0; k < pos.size(); ++k) {
            double rho = std::sqrt(pos[k].x()*pos[k].x() + pos[k].y()*pos[k].y());
            hVz.Fill(pos[k].z());
            hVrho.Fill(rho);
            bool isSec = (st[k] & 0b0110) != 0;  // bits 2 or 3
            bool isDummy = (st[k] & 0b0001) != 0;
            if (!isSec && !isDummy) hVprim.Fill(rho);
            hVno.Fill(nO[k]);
        }
    }

    { TCanvas c("c","",800,600); hNv.Draw(); save(c, out, "09_vertices_per_event"); }
    { TCanvas c("c","",800,600); c.SetLogy(); hVz.Draw(); save(c, out, "10_vertex_z"); }
    { TCanvas c("c","",800,600); c.SetLogy(); hVrho.Draw();
      TLatex txt; txt.SetNDC(); txt.SetTextSize(0.032);
      txt.DrawLatex(0.45, 0.82, "primary + secondaries");
      txt.DrawLatex(0.45, 0.77, "primaries near origin");
      txt.DrawLatex(0.45, 0.72, "secondaries from V^{0} decays");
      txt.DrawLatex(0.45, 0.67, "extend to #rho ~ few cm");
      save(c, out, "11_vertex_rho"); }
    { TCanvas c("c","",800,600); hVprim.Draw(); save(c, out, "12_vertex_rho_primary"); }
    { TCanvas c("c","",800,600); c.SetLogy(); hVno.Draw(); save(c, out, "13_vertex_multiplicity"); }
    { TCanvas c("c","",800,600); c.Divide(1,3);
      c.cd(1); hBx.Draw(); c.cd(2); hBy.Draw(); c.cd(3); hBz.Draw();
      save(c, out, "14_beamspot"); }

    std::cout << "done -- plots under " << out << std::endl;
}
