import ROOT
import numpy as np, dataclasses as dc
import sys
import os
import argparse
import math
import itertools
from array import array
from common_functions import *
from binning_and_selections import *

class MyResponse:

    def __init__(self, reco_tree, gen_tree):
        self._hists = {}
        self._resps = {}
        self._treco = reco_tree
        self._tgen = gen_tree
        self._evt_counter = 0
        self._match_r = 0.05
        self._isALEPH = False

    def setALEPH(self):
        self._isALEPH = True
        self._match_r = 9999

    def writeToFile(self, output):
        fout = ROOT.TFile(output, 'recreate')
        #fout.SetCompressionLevel(9)
        fout.cd()
        for key in self._hists.keys():
            self._hists[key].Write()
        for key in self._resps.keys():
            self._resps[key].Write()

        cov_hist_r = ROOT.TH2D(
            f"covariance_matrix_r", 
            f"EEC r Covariance Matrix",
            self._total_bins, 0.5, self._total_bins + 0.5,
            self._total_bins, 0.5, self._total_bins + 0.5
        )

        cov_hist_z = ROOT.TH2D(
            f"covariance_matrix_z", 
            f"EEC z Covariance Matrix",
            self._total_bins, 0.5, self._total_bins + 0.5,
            self._total_bins, 0.5, self._total_bins + 0.5
        )
        
        for i in range(self._total_bins):
            for j in range(self._total_bins):
                cov_hist_r.SetBinContent(i + 1, j + 1, self._sum_of_eec_products_r[i, j])
                cov_hist_z.SetBinContent(i + 1, j + 1, self._sum_of_eec_products_z[i, j])
        
        # Also save the mean EEC as a histogram
        mean_hist_r = ROOT.TH1D(
            f"mean_eec_r",
            f"Mean EEC r",
            self._total_bins, 0.5, self._total_bins + 0.5
        )

        mean_hist_z = ROOT.TH1D(
            f"mean_eec_z",
            f"Mean EEC z",
            self._total_bins, 0.5, self._total_bins + 0.5
        )
        
        for i in range(self._total_bins):
            mean_hist_r.SetBinContent(i + 1, self._sum_of_eecs_r[i])
            mean_hist_z.SetBinContent(i + 1, self._sum_of_eecs_z[i])

        cov_hist_r.Write()
        mean_hist_r.Write()

        cov_hist_z.Write()
        mean_hist_z.Write()
            
        fout.Close()

    # --------------------------------------------------------------------
    def bookHistograms(self):
    
        # ================================================================
        # 1-D  (name  →  edges)
        # ================================================================
        h1d_defs = {
            "counter"         : np.array([0., 1., 2., 3.], dtype="f"),
    
            # global 1-D spectra
            "reco1d_eec_r"      : rbins,
            "gen1d_eec_r"       : rbins,
            "reco1d_eec_z"      : zbins,
            "gen1d_eec_z"       : zbins,

            "gen_thrust"        : tbins,
            "reco_thrust"       : tbins,
            "gen_thrust_log"    : logtbins,
            "reco_thrust_log"   : logtbins,

            "gen_thrust2"        : tbins2,
            "reco_thrust2"       : tbins2,
            "gen_thrustDelphi"  : tbinsDelphi,
            "reco_thrustDelphi" : tbinsDelphi,
            "gen_thrust_log2"    : logtbins2,
            "reco_thrust_log2"   : logtbins2,

            "reco_e"            : ebins,
            "gen_e"             : ebins,

            "reco_e_match"            : ebins,
            "gen_e_match"             : ebins,
            
            # fake / miss bookkeeping
            "fake_r"          : rbins,
            "fake_z"          : zbins,
            "fake_eij"    : eijbins,
            "miss_r"          : rbins,
            "miss_z"          : zbins,
            "miss_eij"    : eijbins,

            "fake_e"      : ebins,
            "miss_e"      : ebins,

            # By definition, thrust unfolding has no fake or miss
            # as gen- and reco-level dists are calculate with the
            # same event selection and the event selection effects
            # will be corrected in a later step
            # but book the histogram anyway
            "fake_thrust"     : tbins,
            "miss_thrust"     : tbins,
            "fake_thrust_log"     : logtbins,
            "miss_thrust_log"     : logtbins,
            "reco_thrust_theta"   : np.array(np.linspace(0, 180, 181)),
            "gen_thrust_theta"    : np.array(np.linspace(0, 180, 181)),
          }
    
        # ================================================================
        # 2-D  (name → (x_edges , y_edges))
        #  ─── only   E_ij × R   and   E_ij × z
        # ================================================================
        h2d_defs = {
            # ----------  E_ij  ×  R  ------------------------------------
            "reco2d_eij_r"           : (rbins,  eijbins),
            "reco2d_eij_r_match" : (rbins,  eijbins),
            "gen2d_eij_r"            : (rbins,  eijbins),
            "gen2d_eij_r_match"  : (rbins,  eijbins),
    
            # ----------  E_ij  ×  z  ------------------------------------
            "reco2d_eij_z"           : (zbins,  eijbins),
            "reco2d_eij_z_match" : (zbins,  eijbins),
            "gen2d_eij_z"            : (zbins,  eijbins),
            "gen2d_eij_z_match"  : (zbins,  eijbins),
    
            # ---------- 2-D response matrices ---------------------------
            "resp_r"         : (rbins,   rbins),
            "resp_z"         : (zbins,   zbins),
            "resp_eij"   : (eijbins, eijbins),

            "resp_thrust"    : (tbins, tbins),
            "resp_thrust_log": (logtbins, logtbins),

            "resp_thrust2"    : (tbins2, tbins2),
            "resp_thrustDelphi"    : (tbinsDelphi, tbinsDelphi),
            "resp_thrust_log2": (logtbins2, logtbins2),

            "resp_e"         : (ebins, ebins),
        }
    
        # ================================================================
        #  BOOK THEM
        # ================================================================
        self._hists = {}                    
    
        # ---------- 1-D --------------------------------------------------
        self._hists.update({
            name: ROOT.TH1F(name, "", len(edges) - 1, edges)
            for name, edges in h1d_defs.items()
        })
    
        # ---------- 2-D --------------------------------------------------
        self._hists.update({
            name: ROOT.TH2D(name, "", len(xe) - 1, xe,
                            len(ye) - 1, ye)
            for name, (xe, ye) in h2d_defs.items()
        })

        self._template_hist_r=self._hists["reco2d_eij_r"]
        self._template_hist_z=self._hists["reco2d_eij_z"]
        self._nx = self._template_hist_r.GetNbinsX() + 2  # +2 for under/overflow
        self._ny = self._template_hist_r.GetNbinsY() + 2  # +2 for under/overflow  
        self._total_bins = self._nx * self._ny

        self._sum_of_eecs_r=np.zeros(self._total_bins)
        self._sum_of_eec_products_r = np.zeros((self._total_bins, self._total_bins)) 

        self._sum_of_eecs_z=np.zeros(self._total_bins)
        self._sum_of_eec_products_z = np.zeros((self._total_bins, self._total_bins)) 

    def bookResponseMatrices(self):

        respname = 'response2d_eij_r'
        self._resps[respname] = ROOT.RooUnfoldResponse(respname, respname)
        self._resps[respname].Setup(self._hists['reco2d_eij_r'], self._hists['gen2d_eij_r'])

        respname = 'response2d_eij_z'
        self._resps[respname] = ROOT.RooUnfoldResponse(respname, respname)
        self._resps[respname].Setup(self._hists['reco2d_eij_z'], self._hists['gen2d_eij_z'])

        respname = 'response_thrust'
        self._resps[respname] = ROOT.RooUnfoldResponse(respname, respname)
        self._resps[respname].Setup(self._hists['reco_thrust'], self._hists['gen_thrust'])

        respname = 'response_thrust_log'
        self._resps[respname] = ROOT.RooUnfoldResponse(respname, respname)
        self._resps[respname].Setup(self._hists['reco_thrust_log'], self._hists['gen_thrust_log'])

    def bookResponseMatricesNoRooUnfold(self):
        hnd_defs = {
            'response2d_eij_r'   : [ rbins,    eijbins,  rbins,    eijbins ],
            'response2d_eij_z'   : [ zbins,    eijbins,  zbins,    eijbins ],
        }

        for name, edges in hnd_defs.items():
            ndim = len(edges)
            # number of bins per axis
            nbins = array('i', [len(e)-1 for e in edges])

            # build a THnD with arbitrary (uniform) [0,1] ranges:
            h = ROOT.THnD(name, name,
                          ndim,
                          nbins,
                          array('d', [0.0]*ndim),
                          array('d', [1.0]*ndim))

            # now override each axis with the real variable bin edges:
            for ax in range(ndim):
                axis = h.GetAxis(ax)
                axis.Set(len(edges[ax]) - 1,
                         array('d', edges[ax]))  # must be C‐array of doubles

            self._resps[name] = h

        # thrust (1D→1D response) is just a TH2F:
        self._resps["response_thrust"] = ROOT.TH2D(
            "response_thrust",     "",
            len(tbins)-1,    tbins,
            len(tbins)-1,    tbins,
        )
        self._resps["response_thrust_log"] = ROOT.TH2D(
            "response_thrust_log", "",
            len(logtbins)-1, logtbins,
            len(logtbins)-1, logtbins,
        )
        self._resps["response_thrust2"] = ROOT.TH2D(
            "response_thrust2",     "",
            len(tbins2)-1,    tbins2,
            len(tbins2)-1,    tbins2,
        )
        self._resps["response_thrustDelphi"] = ROOT.TH2D(
            "response_thrustDelphi",     "",
            len(tbinsDelphi)-1,    tbinsDelphi,
            len(tbinsDelphi)-1,    tbinsDelphi,
        )
        self._resps["response_thrust_log2"] = ROOT.TH2D(
            "response_thrust_log2", "",
            len(logtbins2)-1, logtbins2,
            len(logtbins2)-1, logtbins2,
        )
    
    def loop(self):

        nevt = self._treco.GetEntries()
        for ievt in range(nevt):

            if ievt % 1000 == 0:
                print(f"Processing event {ievt}/{nevt}")

            treco, tgen = self._treco, self._tgen
            treco.GetEntry(ievt)
            tgen .GetEntry(ievt)

            try:
                E_reco   = treco.Energy
            except:
                E_reco = 91.25

            try:
                E_gen   = tgen.Energy
            except:
                E_gen = 91.25

            get      = lambda *names: (np.asarray(getattr(treco, n)) for n in names)
            get_gen  = lambda *names: (np.asarray(getattr(tgen , n)) for n in names)

            if self._isALEPH:
                px,py,pz,m,q,th,pt,eta,phi,d0,z0,pwflag,hp = get(
                    'px','py','pz','mass','charge','theta','pt','eta','phi','d0','z0','pwflag','highPurity')
        
                px_gen,py_gen,pz_gen,m_gen,q_gen,th_gen,pt_gen,eta_gen,phi_gen,pwflag_gen,hp_gen = get_gen(
                    'px','py','pz','mass','charge','theta','pt','eta','phi','pwflag','highPurity')
            else:
                px,py,pz,m,q,th,pt,eta,phi,d0,z0,pwflag,hp,idx,cspidx = get(
                    'px','py','pz','mass','charge','theta','pt','eta','phi','d0','z0','pwflag','highPurity','index','correspondenceIndex')
        
                px_gen,py_gen,pz_gen,m_gen,q_gen,th_gen,pt_gen,eta_gen,phi_gen,pwflag_gen,hp_gen,idx_gen,cspidx_gen = get_gen(
                    'px','py','pz','mass','charge','theta','pt','eta','phi','pwflag','highPurity','index','correspondenceIndex')

            sphericity = calculate_sphericity_with_fallback(px, py, pz)
            sphericity_gen = calculate_sphericity_with_fallback(px_gen, py_gen, pz_gen)

            # Gen level selection neccessary for DELPHI
            # For ALEPH, this is always true anyway
            e = np.sqrt(px_gen**2 + py_gen**2 + pz_gen**2 + m_gen**2)
            if (np.sum(e) - E_gen) > 0.1: continue

            self._hists['counter'].Fill(0.5)

            if self._isALEPH:
                all_results = apply_track_selection_aleph(px=px, py=py, pz=pz, m=m, q=q, th=th,
                                                          pt=pt, eta=eta, phi=phi, d0=d0, z0=z0, pwflag=pwflag, hp=hp,
                                                          px_gen=px_gen, py_gen=py_gen, pz_gen=pz_gen, m_gen=m_gen, q_gen=q_gen,
                                                          th_gen=th_gen, pt_gen=pt_gen, eta_gen=eta_gen, phi_gen=phi_gen,
                                                          pwflag_gen=pwflag_gen, hp_gen=hp_gen)
                sel_c = all_results['sel_c']          # Reco charged particles
                sel = all_results['sel']              # All reco particles
                sel_c_gen = all_results['sel_c_gen']  # Gen charged particles
                sel_gen = all_results['sel_gen']      # All gen particles
                
                #sel_c, sel, sel_c_gen, sel_gen = apply_track_selection_aleph(
                #    px, py, pz, m, q, th, pt, eta, phi, d0, z0, pwflag, hp,
                #    px_gen, py_gen, pz_gen, m_gen, q_gen, th_gen, pt_gen, eta_gen, phi_gen, pwflag_gen, hp_gen
                #)
            else:
                all_results = apply_track_selection_delphi(px=px, py=py, pz=pz, m=m, q=q, th=th,
                                                           pt=pt, eta=eta, phi=phi, d0=d0, z0=z0, pwflag=pwflag, hp=hp,
                                                           px_gen=px_gen, py_gen=py_gen, pz_gen=pz_gen, m_gen=m_gen, q_gen=q_gen,
                                                           th_gen=th_gen, pt_gen=pt_gen, eta_gen=eta_gen, phi_gen=phi_gen,
                                                           pwflag_gen=pwflag_gen, hp_gen=hp_gen)
                sel_c = all_results['sel_c']          # Reco charged particles
                sel = all_results['sel']              # All reco particles
                sel_c_gen = all_results['sel_c_gen']  # Gen charged particles
                sel_gen = all_results['sel_gen']      # All gen particles
                
                #sel_c, sel, sel_c_gen, sel_gen = apply_track_selection_delphi(
                #    px, py, pz, m, q, th, pt, eta, phi, d0, z0, pwflag, hp,
                #    px_gen, py_gen, pz_gen, m_gen, q_gen, th_gen, pt_gen, eta_gen, phi_gen, pwflag_gen, hp_gen
                #)

            if self._isALEPH:
                px_c, py_c, pz_c, m_c, q_c, pt_c, th_c, phi_c = (
                    v1[sel_c] for v1 in (px, py, pz, m, q, pt, th, phi)
                )

                px_gen_c, py_gen_c, pz_gen_c, m_gen_c, q_gen_c, pt_gen_c, th_gen_c, phi_gen_c = (
                    v2[sel_c_gen] for v2 in (px_gen, py_gen, pz_gen, m_gen, q_gen, pt_gen, th_gen, phi_gen)
                )
            else:
                px_c, py_c, pz_c, m_c, q_c, pt_c, th_c, phi_c, idx_c, cspidx_c = (
                    v1[sel_c] for v1 in (px, py, pz, m, q, pt, th, phi, idx, cspidx)
                )

                px_gen_c, py_gen_c, pz_gen_c, m_gen_c, q_gen_c, pt_gen_c, th_gen_c, phi_gen_c, idx_gen_c, cspidx_gen_c = (
                    v2[sel_c_gen] for v2 in (px_gen, py_gen, pz_gen, m_gen, q_gen, pt_gen, th_gen, phi_gen, idx_gen, cspidx_gen)
                )

            px_n, py_n, pz_n, m_n, q_n, pt_n, th_n, phi_n = (
                v3[sel] for v3 in (px, py, pz, m, q, pt, th, phi)
            )

            px_gen_n, py_gen_n, pz_gen_n, m_gen_n, q_gen_n, pt_gen_n, th_gen_n, phi_gen_n = (
                v4[sel_gen] for v4 in (px_gen, py_gen, pz_gen, m_gen, q_gen, pt_gen, th_gen, phi_gen)
            )

            # --- Event selection and thrust ---
            e_c = np.sqrt(px_c**2 + py_c**2 + pz_c**2 + m_c**2)
            e_n = np.sqrt(px_n**2 + py_n**2 + pz_n**2 + m_n**2)

            e_gen_c = np.sqrt(px_gen_c**2 + py_gen_c**2 + pz_gen_c**2 + m_gen_c**2)
            e_gen_n = np.sqrt(px_gen_n**2 + py_gen_n**2 + pz_gen_n**2 + m_gen_n**2)

            if self._isALEPH:
                rec_c = P4Block.build(px_c ,py_c ,pz_c ,q_c ,
                                        pt_c ,th_c ,phi_c ,e_c ,m_c)
                gen_c = P4Block.build(px_gen_c ,py_gen_c ,pz_gen_c ,q_gen_c ,
                                        pt_gen_c ,th_gen_c ,phi_gen_c ,e_gen_c, m_gen_c)
            else:
                rec_c = P4Block.build(px_c ,py_c ,pz_c ,q_c ,
                                        pt_c ,th_c ,phi_c ,e_c ,m_c ,idx_c, cspidx_c)
                gen_c = P4Block.build(px_gen_c ,py_gen_c ,pz_gen_c ,q_gen_c ,
                                        pt_gen_c ,th_gen_c ,phi_gen_c ,e_gen_c, m_gen_c, idx_gen_c, cspidx_gen_c)               

            rec = P4Block.build(px_n ,py_n ,pz_n ,q_n ,
                                      pt_n ,th_n ,phi_n ,e_n)
            gen = P4Block.build(px_gen_n ,py_gen_n ,pz_gen_n ,q_gen_n ,
                                      pt_gen_n ,th_gen_n ,phi_gen_n ,e_gen_n)

            axis_n_met, T_n_met = thrust_axis_fast(rec.p, include_met=True)
            axis_gen_n_met, T_gen_n_met = thrust_axis_fast(gen.p, include_met=True)

            theta_Tu = thrust_theta(axis_n_met, T_n_met, fold=False)
            theta_gen_Tu = thrust_theta(axis_gen_n_met, T_n_met, fold=False)

            self._hists["reco_thrust_theta"].Fill(np.degrees(theta_Tu))
            self._hists["gen_thrust_theta"].Fill(np.degrees(theta_gen_Tu))

            if self._isALEPH:
                results = apply_event_selection_aleph(
                    e_c, e_n, e_gen_c, e_gen_n, sphericity["cos_theta_v1"], sphericity_gen["cos_theta_v1"]
                )
                pass_reco = results['pass_reco']
                pass_gen = results['pass_reco']
                pass_reco = treco.passesSTheta > 0.5 and treco.passesNTrkMin > 0.5 and treco.passesTotalChgEnergyMin > 0.5 and treco.passesNeuNch > 0.5
            else:
                results = apply_event_selection_delphi(
                    e_c, e_n, e_gen_c, e_gen_n, theta_Tu, theta_gen_Tu, E_reco, E_gen
                )
                pass_reco = results['pass_reco']
                pass_gen = results['pass_reco']

            #if not (pass_reco and pass_gen):
            if not pass_reco:
                continue
            
            self._hists['counter'].Fill(1.5)
            self._evt_counter += 1

            self._hists["reco_thrust"].Fill(T_n_met)
            self._hists["reco_thrust_log"].Fill(np.log(1-T_n_met))
            self._hists["gen_thrust"].Fill(T_gen_n_met)
            self._hists["gen_thrust_log"].Fill(np.log(1-T_gen_n_met))
            self._hists["resp_thrust"].Fill(T_n_met, T_gen_n_met)
            self._hists["resp_thrust_log"].Fill(np.log(1-T_n_met), np.log(1-T_gen_n_met))
            self._resps["response_thrust"].Fill(T_n_met, T_gen_n_met)
            self._resps["response_thrust_log"].Fill(np.log(1-T_n_met), np.log(1-T_gen_n_met))

            self._hists["reco_thrust2"].Fill(T_n_met)
            self._hists["reco_thrustDelphi"].Fill(1-T_n_met)
            self._hists["reco_thrust_log2"].Fill(np.log(1-T_n_met))
            self._hists["gen_thrust2"].Fill(T_gen_n_met)
            self._hists["gen_thrustDelphi"].Fill(1-T_gen_n_met)
            self._hists["gen_thrust_log2"].Fill(np.log(1-T_gen_n_met))
            self._hists["resp_thrust2"].Fill(T_n_met, T_gen_n_met)
            self._hists["resp_thrustDelphi"].Fill(1-T_n_met, 1-T_gen_n_met)
            self._hists["resp_thrust_log2"].Fill(np.log(1-T_n_met), np.log(1-T_gen_n_met))
            self._resps["response_thrust2"].Fill(T_n_met, T_gen_n_met)
            self._resps["response_thrustDelphi"].Fill(1-T_n_met, 1-T_gen_n_met)
            self._resps["response_thrust_log2"].Fill(np.log(1-T_n_met), np.log(1-T_gen_n_met))

            ireco, igen, imiss, ifake = match_angular(rec_c, gen_c, self._match_r)
            #ireco, igen, imiss, ifake = match_correspondence_index(rec_c, gen_c, "trk")

            # ----------------------------------------------------------------------
            # ➊ matched pairs  (i ↔ j both in ireco / igen)
            # ----------------------------------------------------------------------
            for k in range(len(ireco)):
                i_rec = ireco[k]
                i_gen = igen [k]
                self._hists['reco_e'].Fill(rec_c.e[i_rec])
                self._hists['gen_e'].Fill(gen_c.e[i_gen])
                self._hists['reco_e_match'].Fill(rec_c.e[i_rec])
                self._hists['gen_e_match'].Fill(gen_c.e[i_gen])
                self._hists['resp_e'].Fill(rec_c.e[i_rec], gen_c.e[i_gen])
                for l in range(k+1, len(ireco)):          # k < l ⇒ no double count
                    j_rec = ireco[l]
                    j_gen = igen [l]
            
                    # --- observables ------------------------------------------------
                    Eij_rec = rec_c.e[i_rec]*rec_c.e[j_rec] / E_reco**2
                    Eij_gen = gen_c.e[i_gen]*gen_c.e[j_gen] / E_gen **2
            
                    r_rec   = calcAngle(rec_c.p[i_rec],  rec_c.p[j_rec])
                    r_gen   = calcAngle(gen_c.p[i_gen],  gen_c.p[j_gen])
            
                    z_rec   = (1.0 - np.cos(r_rec))/2.0
                    z_gen   = (1.0 - np.cos(r_gen))/2.0
            
                    # --- histograms -------------------------------------------------
                    self._hists['resp_eij']  .Fill(Eij_rec, Eij_gen)
                    self._hists['resp_r']    .Fill(r_rec  , r_gen)
                    self._hists['resp_z']    .Fill(z_rec  , z_gen)
            
                    self._hists['reco2d_eij_r']          .Fill(r_rec, Eij_rec)
                    self._hists['reco2d_eij_r_match'].Fill(r_rec, Eij_rec)
                    self._hists['gen2d_eij_r']           .Fill(r_gen, Eij_gen)
                    self._hists['gen2d_eij_r_match'] .Fill(r_gen, Eij_gen)
            
                    self._hists['reco2d_eij_z']          .Fill(z_rec, Eij_rec)
                    self._hists['reco2d_eij_z_match'].Fill(z_rec, Eij_rec)
                    self._hists['gen2d_eij_z']           .Fill(z_gen, Eij_gen)
                    self._hists['gen2d_eij_z_match'] .Fill(z_gen, Eij_gen)
            
                    self._hists['reco1d_eec_r'].Fill(r_rec, Eij_rec)
                    self._hists['gen1d_eec_r'] .Fill(r_gen, Eij_gen)
                    self._hists['reco1d_eec_z'].Fill(z_rec, Eij_rec)
                    self._hists['gen1d_eec_z'] .Fill(z_gen, Eij_gen)
            
                    # --- responses --------------------------------------------------
                    coords_r = array('d', [r_rec, Eij_rec, r_gen, Eij_gen])
                    self._resps["response2d_eij_r"].Fill(coords_r)
                    coords_z = array('d', [z_rec, Eij_rec, z_gen, Eij_gen])
                    self._resps["response2d_eij_z"].Fill(coords_z)
            
            # ----------------------------------------------------------------------
            # ➋ missed generator pairs (both indices in ‘imiss’ **or** one in imiss)
            # ----------------------------------------------------------------------
            for k, i_gen in enumerate(imiss):
                self._hists['miss_e'].Fill(gen_c.e[i_gen])
                self._hists['gen_e'].Fill(gen_c.e[i_gen])
                for j_gen in imiss[k+1:]:
                    Eij = gen_c.e[i_gen]*gen_c.e[j_gen] / E_gen**2
                    r   = calcAngle(gen_c.p[i_gen], gen_c.p[j_gen])
                    z   = (1.0 - np.cos(r))/2.0
            
                    self._hists['miss_eij'] .Fill(Eij)
                    self._hists['miss_r']   .Fill(r, Eij)
                    self._hists['miss_z']   .Fill(z, Eij)
            
                    self._hists['gen2d_eij_r'].Fill(r, Eij)
                    self._hists['gen2d_eij_z'].Fill(z, Eij)
                    self._hists['gen1d_eec_r'].Fill(r, Eij)
                    self._hists['gen1d_eec_z'].Fill(z, Eij)
            
                # mix : one index missed, one matched
                for j_gen in igen:
                    if i_gen == j_gen:          # safety
                        continue
                    Eij = gen_c.e[i_gen]*gen_c.e[j_gen] / E_gen**2
                    r   = calcAngle(gen_c.p[i_gen], gen_c.p[j_gen])
                    z   = (1.0 - np.cos(r))/2.0
            
                    self._hists['miss_eij'] .Fill(Eij)
                    self._hists['miss_r']   .Fill(r, Eij)
                    self._hists['miss_z']   .Fill(z, Eij)
            
                    self._hists['gen2d_eij_r'].Fill(r, Eij)
                    self._hists['gen2d_eij_z'].Fill(z, Eij)
                    self._hists['gen1d_eec_r'].Fill(r, Eij)
                    self._hists['gen1d_eec_z'].Fill(z, Eij)
            
            # ----------------------------------------------------------------------
            # ➌ fake reconstructed pairs (both in ‘ifake’ **or** one in ifake)
            # ----------------------------------------------------------------------
            event_pairs_r = []
            event_pairs_z = []
            for k, i_rec in enumerate(ifake):
                self._hists['fake_e'].Fill(rec_c.e[i_rec])
                self._hists['reco_e'].Fill(rec_c.e[i_rec])
                for j_rec in ifake[k+1:]:
                    Eij = rec_c.e[i_rec]*rec_c.e[j_rec] / E_reco**2
                    r   = calcAngle(rec_c.p[i_rec], rec_c.p[j_rec])
                    z   = (1.0 - np.cos(r))/2.0
            
                    self._hists['fake_eij'].Fill(Eij)
                    self._hists['fake_r']  .Fill(r, Eij)
                    self._hists['fake_z']  .Fill(z, Eij)
            
                    self._hists['reco2d_eij_r'].Fill(r, Eij)
                    self._hists['reco2d_eij_z'].Fill(z, Eij)
                    self._hists['reco1d_eec_r'].Fill(r, Eij)
                    self._hists['reco1d_eec_z'].Fill(z, Eij)

                    event_pairs_r.append((r, Eij, 1))
                    event_pairs_z.append((z, Eij, 1))
            
                # mix : one index fake, one matched
                for j_rec in ireco:
                    if i_rec == j_rec:
                        continue
                    Eij = rec_c.e[i_rec]*rec_c.e[j_rec] / E_reco**2
                    r   = calcAngle(rec_c.p[i_rec], rec_c.p[j_rec])
                    z   = (1.0 - np.cos(r))/2.0
            
                    self._hists['fake_eij'].Fill(Eij)
                    self._hists['fake_r']  .Fill(r, Eij)
                    self._hists['fake_z']  .Fill(z, Eij)
                    
                    self._hists['reco1d_eec_r'].Fill(r, Eij)
                    self._hists['reco1d_eec_z'].Fill(z, Eij)
                    self._hists['reco2d_eij_r'].Fill(r, Eij)
                    self._hists['reco2d_eij_z'].Fill(z, Eij)

                    event_pairs_r.append((r, Eij, 1))
                    event_pairs_z.append((z, Eij, 1))

            event_eec_r = calculate_event_eec_histogram(event_pairs_r, self._template_hist_r, self._nx, self._ny)
            event_eec_z = calculate_event_eec_histogram(event_pairs_z, self._template_hist_z, self._nx, self._ny)

            self._sum_of_eecs_r += event_eec_r
            self._sum_of_eec_products_r += np.outer(event_eec_r, event_eec_r) 

            self._sum_of_eecs_z += event_eec_z
            self._sum_of_eec_products_z += np.outer(event_eec_z, event_eec_z)        

    def findBin(self, bin_edges, value):
        bin_index = np.digitize(value, bin_edges) - 1

        return bin_index

    def normalize(self):
        self._hists['reco1d_eec_r'].Scale(1./self._evt_counter)
        self._hists['reco1d_eec_r']=normalizeByBinWidth(self._hists['reco1d_eec_r'])
        
        self._hists['gen1d_eec_r'].Scale(1./self._evt_counter)
        self._hists['gen1d_eec_r']=normalizeByBinWidth(self._hists['gen1d_eec_r'])

        self._hists['reco1d_eec_z'].Scale(1./self._evt_counter)
        self._hists['reco1d_eec_z']=normalizeByBinWidth(self._hists['reco1d_eec_z'])
        
        self._hists['gen1d_eec_z'].Scale(1./self._evt_counter)
        self._hists['gen1d_eec_z']=normalizeByBinWidth(self._hists['gen1d_eec_z'])
    

if __name__ == "__main__":

    #filename = '/eos/user/z/zhangj/DELPHI/simulation/v94c/91.25/kk2f4146_qqpy/nanoaod_kk2f4146_qqpy_91.25_40001.sdst.root'
    filename = '/eos/user/z/zhangj/ALEPH/SamplesLEP1/ALEPHMC/LEP1MC1994_recons_aftercut-020.root'
    filenameout = "response_test.root"

    parser = argparse.ArgumentParser()
    parser.add_argument("infiles", nargs='?', default=filename, help="name of input files")
    parser.add_argument("outfile", nargs='?', default=filenameout, help="name of input files")
    args = parser.parse_args()

    treco = 't'
    if "ALEPH" in args.infiles:
        tgen='tgen'
    else:
        tgen = 'tgenBefore'

    t_reco = ROOT.TChain(treco)
    t_gen = ROOT.TChain(tgen)

    print("Reading input from: ",args.infiles)
    InputRootFiles=[]
    if args.infiles.find(".root")>-1:
        InputRootFiles.append(args.infiles)
    else:
        ## read from list
        InputRootFiles=ReadFilesFromList(args.infiles)
    
    for f in InputRootFiles:
        t_reco.Add(f)
        t_gen.Add(f)

    fnameout = args.outfile

    response = MyResponse(t_reco, t_gen)
    if "ALEPH" in args.infiles:
        response.setALEPH()
        print("ALEPH")
    response.bookHistograms()
    response.bookResponseMatricesNoRooUnfold()
    response.loop()
    response.writeToFile(fnameout)