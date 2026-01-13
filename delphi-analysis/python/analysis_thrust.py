#!/usr/bin/env python3
"""
Optimized Thrust Analysis
-------------------------
Single-pass thrust analysis with systematic variations.

For MC:
  - Fills gen-level histograms (before selection, no systematics)
  - Fills reco-level histograms (after selection AND event selection, all systematics)
  - Fills response matrices (after selection AND event selection, all systematics)

For Data:
  - Fills reco-level histograms (after selection AND event selection, all systematics)

Usage:
  MC:   python analysis_thrust.py input.root output.root --is-mc
  Data: python analysis_thrust.py input.root output.root
  Single systematic: python analysis_thrust.py input.root output.root --systematic nominal
"""

import ROOT
import numpy as np
import argparse
import math
from array import array

# Import from common modules
from common_functions import (
    thrust_axis_fast_optimized,
    thrust_theta,
    apply_momentum_bias,
    apply_fake_drop_charged,
    randomly_drop_particles,
    calc_multiplicity_weight_linear,
    missing_p,
)
from binning_and_selections import (
    tbins2,
    logtbins2,
    apply_track_selection_delphi,
    apply_event_selection_delphi,
)

# ============================================================================
# SYSTEMATIC CONFIGURATION
# ============================================================================

SYSTEMATICS = {
    'nominal': {
        'type': 'selection',
        'charged_pt_min': 0.4,
        'neutral_e_min': 0.5,
    },
    
    # Momentum scale systematics
    'charged_scale_up': {
        'type': 'momentum',
        'momentum_scale_charged': 1.01,
    },
    'charged_scale_down': {
        'type': 'momentum',
        'momentum_scale_charged': 0.99,
    },
    'neutral_scale_up': {
        'type': 'momentum',
        'momentum_scale_neutral': 1.01,
    },
    'neutral_scale_down': {
        'type': 'momentum',
        'momentum_scale_neutral': 0.99,
    },
    
    # Efficiency systematics
    'fake_drop': {
        'type': 'efficiency',
        'apply_fake_drop': True,
    },
    'charged_eff': {
        'type': 'efficiency',
        'drop_charged_fraction': 0.02,
    },
    'neutral_eff': {
        'type': 'efficiency',
        'drop_neutral_fraction': 0.02,
    },
    
    # Event multiplicity reweighting
    'evt_weight': {
        'type': 'reweighting',
        'use_evt_weight': True,
    },
    
    # Selection systematics
    'neutral_e_up': {
        'type': 'selection',
        'charged_pt_min': 0.4,
        'neutral_e_min': 5.0,
    },
    'neutral_e_down': {
        'type': 'selection',
        'charged_pt_min': 0.4,
        'neutral_e_min': 0.3,
    },
    'charged_pt_up': {
        'type': 'selection',
        'charged_pt_min': 0.5,
        'neutral_e_min': 0.5,
    },
    'charged_pt_down': {
        'type': 'selection',
        'charged_pt_min': 0.3,
        'neutral_e_min': 0.5,
    },
}

# ============================================================================
# HISTOGRAM CREATION FUNCTIONS
# ============================================================================

def create_gen_histograms():
    """
    Create gen-level histograms (before selection, no systematics).
    Only used for MC.
    """
    hists = {}
    
    # Standard binning (tbins2)
    hists['Thrust_before2'] = ROOT.TH1D('Thrust_before2', '', len(tbins2)-1, array('d', tbins2))
    hists['Thrust_before2_Escheme'] = ROOT.TH1D('Thrust_before2_Escheme', '', len(tbins2)-1, array('d', tbins2))
    hists['ThrustC_before2'] = ROOT.TH1D('ThrustC_before2', '', len(tbins2)-1, array('d', tbins2))
    hists['ThrustC_before2_Escheme'] = ROOT.TH1D('ThrustC_before2_Escheme', '', len(tbins2)-1, array('d', tbins2))
    
    # Log binning (logtbins2)
    hists['Thrust_before_log2'] = ROOT.TH1D('Thrust_before_log2', '', len(logtbins2)-1, array('d', logtbins2))
    hists['Thrust_before_log2_Escheme'] = ROOT.TH1D('Thrust_before_log2_Escheme', '', len(logtbins2)-1, array('d', logtbins2))
    hists['ThrustC_before_log2'] = ROOT.TH1D('ThrustC_before_log2', '', len(logtbins2)-1, array('d', logtbins2))
    hists['ThrustC_before_log2_Escheme'] = ROOT.TH1D('ThrustC_before_log2_Escheme', '', len(logtbins2)-1, array('d', logtbins2))
    
    return hists


def create_reco_histograms(syst_name):
    """
    Create reco-level histograms for one systematic.
    Used for both MC and data.
    """
    hists = {}
    
    # Standard binning (tbins2)
    hists[f'ThrustMissPNC2_{syst_name}'] = ROOT.TH1D(f'ThrustMissPNC2_{syst_name}', '', 
                                                      len(tbins2)-1, array('d', tbins2))
    hists[f'ThrustMissPNC2_Escheme_{syst_name}'] = ROOT.TH1D(f'ThrustMissPNC2_Escheme_{syst_name}', '', 
                                                              len(tbins2)-1, array('d', tbins2))
    hists[f'ThrustC2_{syst_name}'] = ROOT.TH1D(f'ThrustC2_{syst_name}', '', 
                                                len(tbins2)-1, array('d', tbins2))
    hists[f'ThrustC2_Escheme_{syst_name}'] = ROOT.TH1D(f'ThrustC2_Escheme_{syst_name}', '', 
                                                        len(tbins2)-1, array('d', tbins2))
    
    # Log binning (logtbins2)
    hists[f'ThrustMissPNCLog2_{syst_name}'] = ROOT.TH1D(f'ThrustMissPNCLog2_{syst_name}', '', 
                                                         len(logtbins2)-1, array('d', logtbins2))
    hists[f'ThrustMissPNCLog2_Escheme_{syst_name}'] = ROOT.TH1D(f'ThrustMissPNCLog2_Escheme_{syst_name}', '', 
                                                                 len(logtbins2)-1, array('d', logtbins2))
    hists[f'ThrustCLog2_{syst_name}'] = ROOT.TH1D(f'ThrustCLog2_{syst_name}', '', 
                                                   len(logtbins2)-1, array('d', logtbins2))
    hists[f'ThrustCLog2_Escheme_{syst_name}'] = ROOT.TH1D(f'ThrustCLog2_Escheme_{syst_name}', '', 
                                                           len(logtbins2)-1, array('d', logtbins2))
    
    return hists


def create_response_matrices(syst_name):
    """
    Create 2D response matrices for one systematic.
    For MC only - all systematics get response matrices.
    """
    matrices = {}
    
    # Standard binning (tbins2 x tbins2)
    matrices[f'response_thrust2_{syst_name}'] = ROOT.TH2D(
        f'response_thrust2_{syst_name}', '',
        len(tbins2)-1, array('d', tbins2),
        len(tbins2)-1, array('d', tbins2)
    )
    matrices[f'response_thrust2_Escheme_{syst_name}'] = ROOT.TH2D(
        f'response_thrust2_Escheme_{syst_name}', '',
        len(tbins2)-1, array('d', tbins2),
        len(tbins2)-1, array('d', tbins2)
    )
    matrices[f'response_thrustC2_{syst_name}'] = ROOT.TH2D(
        f'response_thrustC2_{syst_name}', '',
        len(tbins2)-1, array('d', tbins2),
        len(tbins2)-1, array('d', tbins2)
    )
    matrices[f'response_thrustC2_Escheme_{syst_name}'] = ROOT.TH2D(
        f'response_thrustC2_Escheme_{syst_name}', '',
        len(tbins2)-1, array('d', tbins2),
        len(tbins2)-1, array('d', tbins2)
    )
    
    # Log binning (logtbins2 x logtbins2)
    matrices[f'response_thrust_log2_{syst_name}'] = ROOT.TH2D(
        f'response_thrust_log2_{syst_name}', '',
        len(logtbins2)-1, array('d', logtbins2),
        len(logtbins2)-1, array('d', logtbins2)
    )
    matrices[f'response_thrust_log2_Escheme_{syst_name}'] = ROOT.TH2D(
        f'response_thrust_log2_Escheme_{syst_name}', '',
        len(logtbins2)-1, array('d', logtbins2),
        len(logtbins2)-1, array('d', logtbins2)
    )
    matrices[f'response_thrustC_log2_{syst_name}'] = ROOT.TH2D(
        f'response_thrustC_log2_{syst_name}', '',
        len(logtbins2)-1, array('d', logtbins2),
        len(logtbins2)-1, array('d', logtbins2)
    )
    matrices[f'response_thrustC_log2_Escheme_{syst_name}'] = ROOT.TH2D(
        f'response_thrustC_log2_Escheme_{syst_name}', '',
        len(logtbins2)-1, array('d', logtbins2),
        len(logtbins2)-1, array('d', logtbins2)
    )
    
    return matrices

# ============================================================================
# SYSTEMATIC VARIATION FUNCTIONS
# ============================================================================

def apply_systematic_variation(px, py, pz, m, q, pt, th, d0, z0, config):
    """
    Apply systematic variations to particle momenta and efficiency.
    Each systematic varies only ONE thing.
    
    Parameters:
    -----------
    px, py, pz, m, q, pt, th, d0, z0 : arrays
        Particle properties
    config : dict
        Systematic configuration
    
    Returns:
    --------
    dict with all particle arrays (dropped particles removed)
    """
    # Make copies
    px_var = px.copy()
    py_var = py.copy()
    pz_var = pz.copy()
    m_var = m.copy()
    q_var = q.copy()
    pt_var = pt.copy()
    th_var = th.copy()
    d0_var = d0.copy()
    z0_var = z0.copy()
    
    # Identify charged and neutral particles
    is_charged = np.abs(q_var) > 0.1
    
    # === Momentum scale systematics ===
    if config.get('momentum_scale_charged'):
        scale = config['momentum_scale_charged']
        px_var[is_charged] *= scale
        py_var[is_charged] *= scale
        pz_var[is_charged] *= scale
        # Update pt and theta
        pt_var[is_charged] = np.sqrt(px_var[is_charged]**2 + py_var[is_charged]**2)
        th_var[is_charged] = np.arctan2(pt_var[is_charged], pz_var[is_charged])
    
    if config.get('momentum_scale_neutral'):
        scale = config['momentum_scale_neutral']
        px_var[~is_charged] *= scale
        py_var[~is_charged] *= scale
        pz_var[~is_charged] *= scale
        # Update theta for neutrals
        pt_neutral = np.sqrt(px_var[~is_charged]**2 + py_var[~is_charged]**2)
        th_var[~is_charged] = np.arctan2(pt_neutral, pz_var[~is_charged])
    
    # === Efficiency systematics ===
    keep_mask = np.ones(len(px_var), dtype=bool)
    
    if config.get('apply_fake_drop'):
        charged_indices = np.where(is_charged)[0]
        if len(charged_indices) > 0:
            fake_dropped = apply_fake_drop_charged(
                px_var[is_charged], 
                py_var[is_charged], 
                pz_var[is_charged], 
                m_var[is_charged],
                seed=config.get('fake_drop_seed')
            )
            keep_mask[charged_indices] = fake_dropped['mask']
    
    if config.get('drop_charged_fraction', 0.0) > 0:
        charged_indices = np.where(is_charged)[0]
        if len(charged_indices) > 0:
            dropped = randomly_drop_particles(
                px_var[is_charged],
                py_var[is_charged],
                pz_var[is_charged],
                m_var[is_charged],
                drop_fraction=config['drop_charged_fraction'],
                seed=config.get('drop_charged_seed')
            )
            keep_mask[charged_indices] = dropped['mask']
    
    if config.get('drop_neutral_fraction', 0.0) > 0:
        neutral_indices = np.where(~is_charged)[0]
        if len(neutral_indices) > 0:
            dropped = randomly_drop_particles(
                px_var[~is_charged],
                py_var[~is_charged],
                pz_var[~is_charged],
                m_var[~is_charged],
                drop_fraction=config['drop_neutral_fraction'],
                seed=config.get('drop_neutral_seed')
            )
            keep_mask[neutral_indices] = dropped['mask']
    
    # Apply keep mask to all arrays
    return {
        'px': px_var[keep_mask],
        'py': py_var[keep_mask],
        'pz': pz_var[keep_mask],
        'm': m_var[keep_mask],
        'q': q_var[keep_mask],
        'pt': pt_var[keep_mask],
        'th': th_var[keep_mask],
        'd0': d0_var[keep_mask],
        'z0': z0_var[keep_mask],
    }


def calculate_all_thrust_variants(px, py, pz, m, q, sel_c, sel_n, include_met=True):
    """
    Calculate all thrust variants for one event.
    
    Parameters:
    -----------
    px, py, pz, m, q : arrays
        Particle properties (already selected)
    sel_c, sel_n : boolean arrays
        Masks for charged and neutral particles
    include_met : bool
        Whether to include missing momentum (True for reco, False for gen)
    
    Returns:
    --------
    dict with thrust values, log values, axis vectors, and T values
    """
    # Recalculate energy from momentum
    e = np.sqrt(px**2 + py**2 + pz**2 + m**2)
    
    # === Standard (P-scheme) Thrust ===
    p3_all = np.stack([px, py, pz], axis=1)
    axis_all, T_all = thrust_axis_fast_optimized(p3_all, include_met=include_met)
    thrust_val = 1 - T_all
    thrust_log_val = np.log(thrust_val) if thrust_val > 0 else -10.0
    
    # === E-scheme Thrust ===
    p3_norm = np.linalg.norm(p3_all, axis=1, keepdims=True)
    n = p3_all / np.where(p3_norm > 0, p3_norm, 1.0)
    p3_escheme = e[:, np.newaxis] * n
    axis_escheme, T_escheme = thrust_axis_fast_optimized(p3_escheme, include_met=include_met)
    thrust_escheme_val = 1 - T_escheme
    thrust_log_escheme_val = np.log(thrust_escheme_val) if thrust_escheme_val > 0 else -10.0
    
    # === Charged-only Thrust ===
    if np.sum(sel_c) > 0:
        p3_charged = np.stack([px[sel_c], py[sel_c], pz[sel_c]], axis=1)
        axis_charged, T_charged = thrust_axis_fast_optimized(p3_charged, include_met=include_met)
        thrust_charged_val = 1 - T_charged
        thrust_charged_log_val = np.log(thrust_charged_val) if thrust_charged_val > 0 else -10.0
    else:
        axis_charged = np.zeros(3)
        T_charged = 0.0
        thrust_charged_val = 0.0
        thrust_charged_log_val = -10.0
    
    # === Charged-only E-scheme Thrust ===
    if np.sum(sel_c) > 0:
        e_charged = e[sel_c]
        p3_charged_norm = np.linalg.norm(p3_charged, axis=1, keepdims=True)
        n_charged = p3_charged / np.where(p3_charged_norm > 0, p3_charged_norm, 1.0)
        p3_charged_escheme = e_charged[:, np.newaxis] * n_charged
        axis_charged_escheme, T_charged_escheme = thrust_axis_fast_optimized(
            p3_charged_escheme, include_met=include_met
        )
        thrust_charged_escheme_val = 1 - T_charged_escheme
        thrust_charged_log_escheme_val = np.log(thrust_charged_escheme_val) if thrust_charged_escheme_val > 0 else -10.0
    else:
        axis_charged_escheme = np.zeros(3)
        T_charged_escheme = 0.0
        thrust_charged_escheme_val = 0.0
        thrust_charged_log_escheme_val = -10.0
    
    return {
        # Thrust values
        'thrust': thrust_val,
        'thrust_log': thrust_log_val,
        'thrust_escheme': thrust_escheme_val,
        'thrust_log_escheme': thrust_log_escheme_val,
        'thrust_charged': thrust_charged_val,
        'thrust_charged_log': thrust_charged_log_val,
        'thrust_charged_escheme': thrust_charged_escheme_val,
        'thrust_charged_log_escheme': thrust_charged_log_escheme_val,
        # Axis and T values for event selection
        'axis': axis_all,
        'T': T_all,
        'axis_escheme': axis_escheme,
        'T_escheme': T_escheme,
        'axis_charged': axis_charged,
        'T_charged': T_charged,
        'axis_charged_escheme': axis_charged_escheme,
        'T_charged_escheme': T_charged_escheme,
    }

# ============================================================================
# EVENT PROCESSING FUNCTIONS
# ============================================================================

def process_gen_level(tree_gen_before, tree_gen, event_idx, gen_hists):
    """
    Process generator-level particles (MC only).
    No event selection applied at gen level.
    
    Parameters:
    -----------
    tree_gen_before : ROOT TTree
        Primary gen tree (tgenBefore)
    tree_gen : ROOT TTree
        Fallback gen tree (tgen) 
    event_idx : int
        Event index
    gen_hists : dict
        Dictionary of gen-level histograms
    
    Returns:
    --------
    dict with gen thrust values (for response matrix filling), or None if event invalid
    """
    # Load from primary tree
    tree_gen_before.GetEntry(event_idx)
    
    # Energy cut
    E = tree_gen_before.Energy
    if abs(E - 91.25) > 1:
        return None
    
    # Load particles
    px = np.array(tree_gen_before.px)
    py = np.array(tree_gen_before.py)
    pz = np.array(tree_gen_before.pz)
    m = np.array(tree_gen_before.mass)
    q = np.array(tree_gen_before.charge)
    
    if len(px) == 0:
        return None
    
    # Check energy conservation
    e = np.sqrt(px**2 + py**2 + pz**2 + m**2)
    
    if abs(np.sum(e) - E) > 0.1:
        # Energy not conserved - switch to fallback tree for THIS event
        tree_gen.GetEntry(event_idx)
        
        # Reload ALL particle data from fallback tree
        px = np.array(tree_gen.px)
        py = np.array(tree_gen.py)
        pz = np.array(tree_gen.pz)
        m = np.array(tree_gen.mass)
        q = np.array(tree_gen.charge)
        
        # Recalculate energy with new data
        e = np.sqrt(px**2 + py**2 + pz**2 + m**2)
    
    if len(px) == 0:
        return None
    
    # Identify charged and neutral (no selection applied for gen)
    is_charged = np.abs(q) > 0.1
    sel_c_gen = is_charged
    sel_n_gen = ~is_charged
    
    # Calculate all thrust variants (no MET for gen level)
    thrust_vals = calculate_all_thrust_variants(
        px, py, pz, m, q,
        sel_c_gen, sel_n_gen,
        include_met=False
    )
    
    # Fill gen histograms (no event selection at gen level)
    gen_hists['Thrust_before2'].Fill(thrust_vals['thrust'])
    gen_hists['Thrust_before_log2'].Fill(thrust_vals['thrust_log'])
    gen_hists['Thrust_before2_Escheme'].Fill(thrust_vals['thrust_escheme'])
    gen_hists['Thrust_before_log2_Escheme'].Fill(thrust_vals['thrust_log_escheme'])
    
    gen_hists['ThrustC_before2'].Fill(thrust_vals['thrust_charged'])
    gen_hists['ThrustC_before_log2'].Fill(thrust_vals['thrust_charged_log'])
    gen_hists['ThrustC_before2_Escheme'].Fill(thrust_vals['thrust_charged_escheme'])
    gen_hists['ThrustC_before_log2_Escheme'].Fill(thrust_vals['thrust_charged_log_escheme'])
    
    return thrust_vals


def process_reco_level(tree_reco, event_idx, systematics, reco_hists, response_matrices, 
                       thrust_gen=None, is_mc=False):
    """
    Process reco-level particles for all systematics.
    Applies event selection separately for each thrust variant.
    
    Parameters:
    -----------
    tree_reco : ROOT TTree
        Reco tree (t)
    event_idx : int
        Event index
    systematics : dict
        Dictionary of systematic configurations to process
    reco_hists : dict
        Dictionary of reco histograms (keyed by systematic name)
    response_matrices : dict
        Dictionary of response matrices (keyed by systematic name)
    thrust_gen : dict, optional
        Gen-level thrust values (for response matrix filling)
    is_mc : bool
        Whether this is MC (affects response matrix filling)
    """
    # Load event
    tree_reco.GetEntry(event_idx)
    
    # Energy cut
    E = tree_reco.Energy
    if abs(E - 91.25) > 1:
        return
    
    # Load reco particles (only essential variables)
    px = np.array(tree_reco.px)
    py = np.array(tree_reco.py)
    pz = np.array(tree_reco.pz)
    m = np.array(tree_reco.mass)
    q = np.array(tree_reco.charge)
    pt = np.array(tree_reco.pt)
    th = np.array(tree_reco.theta)
    d0 = np.array(tree_reco.d0)
    z0 = np.array(tree_reco.z0)
    
    if len(px) == 0:
        return
    
    # Loop over all systematics
    for syst_name, syst_config in systematics.items():
        
        # Apply systematic variations (momentum scale, efficiency drops)
        varied = apply_systematic_variation(
            px, py, pz, m, q, pt, th, d0, z0,
            syst_config
        )
        
        px_var = varied['px']
        py_var = varied['py']
        pz_var = varied['pz']
        m_var = varied['m']
        q_var = varied['q']
        pt_var = varied['pt']
        th_var = varied['th']
        d0_var = varied['d0']
        z0_var = varied['z0']
        
        if len(px_var) == 0:
            continue
        
        # Apply particle selection with configurable cuts
        selection = apply_track_selection_delphi(
            px=px_var, py=py_var, pz=pz_var, m=m_var, q=q_var,
            th=th_var, pt=pt_var,
            d0=d0_var, z0=z0_var,
            charged_pt_min=syst_config.get('charged_pt_min', 0.4),
            neutral_e_min=syst_config.get('neutral_e_min', 0.5),
        )
        
        sel_c = selection['sel_c']
        sel_n = selection['sel_n']
        sel_all = selection['sel']
        
        if not np.any(sel_all):
            continue
        
        # Extract selected particles
        px_sel = px_var[sel_all]
        py_sel = py_var[sel_all]
        pz_sel = pz_var[sel_all]
        m_sel = m_var[sel_all]
        q_sel = q_var[sel_all]
        
        sel_c_final = sel_c[sel_all]
        sel_n_final = sel_n[sel_all]
        
        # Calculate event weight
        evt_weight = 1.0
        if syst_config.get('use_evt_weight'):
            n_charged_sel = np.sum(sel_c_final)
            evt_weight = calc_multiplicity_weight_linear(n_charged_sel)
        
        # Calculate thrust
        thrust_vals = calculate_all_thrust_variants(
            px_sel, py_sel, pz_sel, m_sel, q_sel,
            sel_c_final, sel_n_final,
            include_met=True
        )
        
        # Prepare energies for event selection
        e_c = np.sqrt(px_sel[sel_c_final]**2 + py_sel[sel_c_final]**2 + 
                      pz_sel[sel_c_final]**2 + m_sel[sel_c_final]**2)
        e_all = np.sqrt(px_sel**2 + py_sel**2 + pz_sel**2 + m_sel**2)
        
        # ===================================================================
        # P-SCHEME (ALL PARTICLES) EVENT SELECTION
        # ===================================================================
        theta_thrust_p = thrust_theta(thrust_vals['axis'], thrust_vals['T'], fold=False)
        event_sel_p = apply_event_selection_delphi(
            e_c=e_c, 
            e_n=e_all,
            theta_Tu=theta_thrust_p,
            E_reco=E
        )
        
        if event_sel_p['pass_reco']:
            # Fill P-scheme histograms
            reco_hists[syst_name][f'ThrustMissPNC2_{syst_name}'].Fill(thrust_vals['thrust'], evt_weight)
            reco_hists[syst_name][f'ThrustMissPNCLog2_{syst_name}'].Fill(thrust_vals['thrust_log'], evt_weight)
            
            # Fill P-scheme response matrices
            if is_mc and thrust_gen is not None:
                matrices = response_matrices[syst_name]
                matrices[f'response_thrust2_{syst_name}'].Fill(thrust_vals['thrust'], thrust_gen['thrust'])
                matrices[f'response_thrust_log2_{syst_name}'].Fill(thrust_vals['thrust_log'], thrust_gen['thrust_log'])
        
        # ===================================================================
        # E-SCHEME (ALL PARTICLES) EVENT SELECTION
        # ===================================================================
        theta_thrust_e = thrust_theta(thrust_vals['axis_escheme'], thrust_vals['T_escheme'], fold=False)
        event_sel_e = apply_event_selection_delphi(
            e_c=e_c, 
            e_n=e_all,
            theta_Tu=theta_thrust_e,
            E_reco=E
        )
        
        if event_sel_e['pass_reco']:
            # Fill E-scheme histograms
            reco_hists[syst_name][f'ThrustMissPNC2_Escheme_{syst_name}'].Fill(thrust_vals['thrust_escheme'], evt_weight)
            reco_hists[syst_name][f'ThrustMissPNCLog2_Escheme_{syst_name}'].Fill(thrust_vals['thrust_log_escheme'], evt_weight)
            
            # Fill E-scheme response matrices
            if is_mc and thrust_gen is not None:
                matrices = response_matrices[syst_name]
                matrices[f'response_thrust2_Escheme_{syst_name}'].Fill(thrust_vals['thrust_escheme'], thrust_gen['thrust_escheme'])
                matrices[f'response_thrust_log2_Escheme_{syst_name}'].Fill(thrust_vals['thrust_log_escheme'], thrust_gen['thrust_log_escheme'])
        
        # ===================================================================
        # CHARGED P-SCHEME EVENT SELECTION
        # ===================================================================
        theta_thrust_c = thrust_theta(thrust_vals['axis_charged'], thrust_vals['T_charged'], fold=False)
        event_sel_c = apply_event_selection_delphi(
            e_c=e_c, 
            e_n=e_all,
            theta_Tu=theta_thrust_c,
            E_reco=E
        )
        
        if event_sel_c['pass_reco']:
            # Fill charged P-scheme histograms
            reco_hists[syst_name][f'ThrustC2_{syst_name}'].Fill(thrust_vals['thrust_charged'], evt_weight)
            reco_hists[syst_name][f'ThrustCLog2_{syst_name}'].Fill(thrust_vals['thrust_charged_log'], evt_weight)
            
            # Fill charged P-scheme response matrices
            if is_mc and thrust_gen is not None:
                matrices = response_matrices[syst_name]
                matrices[f'response_thrustC2_{syst_name}'].Fill(thrust_vals['thrust_charged'], thrust_gen['thrust_charged'])
                matrices[f'response_thrustC_log2_{syst_name}'].Fill(thrust_vals['thrust_charged_log'], thrust_gen['thrust_charged_log'])
        
        # ===================================================================
        # CHARGED E-SCHEME EVENT SELECTION
        # ===================================================================
        theta_thrust_ce = thrust_theta(thrust_vals['axis_charged_escheme'], thrust_vals['T_charged_escheme'], fold=False)
        event_sel_ce = apply_event_selection_delphi(
            e_c=e_c, 
            e_n=e_all,
            theta_Tu=theta_thrust_ce,
            E_reco=E
        )
        
        if event_sel_ce['pass_reco']:
            # Fill charged E-scheme histograms
            reco_hists[syst_name][f'ThrustC2_Escheme_{syst_name}'].Fill(thrust_vals['thrust_charged_escheme'], evt_weight)
            reco_hists[syst_name][f'ThrustCLog2_Escheme_{syst_name}'].Fill(thrust_vals['thrust_charged_log_escheme'], evt_weight)
            
            # Fill charged E-scheme response matrices
            if is_mc and thrust_gen is not None:
                matrices = response_matrices[syst_name]
                matrices[f'response_thrustC2_Escheme_{syst_name}'].Fill(thrust_vals['thrust_charged_escheme'], thrust_gen['thrust_charged_escheme'])
                matrices[f'response_thrustC_log2_Escheme_{syst_name}'].Fill(thrust_vals['thrust_charged_log_escheme'], thrust_gen['thrust_charged_log_escheme'])

# ============================================================================
# MAIN FUNCTION
# ============================================================================

def main():
    """
    Main function to run optimized thrust analysis.
    """
    # Parse command line arguments
    parser = argparse.ArgumentParser(
        description='Optimized single-pass thrust analysis with systematic variations'
    )
    parser.add_argument('infile', help='Input ROOT file')
    parser.add_argument('outfile', help='Output ROOT file')
    parser.add_argument('--is-mc', action='store_true', 
                        help='Process as MC (fills gen histograms and response matrices)')
    parser.add_argument('--systematic', default='all',
                        help='Run specific systematic or "all" (default: all)')
    
    args = parser.parse_args()
    
    # Determine which systematics to run
    if args.systematic == 'all':
        systematics_to_run = SYSTEMATICS
    elif args.systematic in SYSTEMATICS:
        systematics_to_run = {args.systematic: SYSTEMATICS[args.systematic]}
    else:
        print(f"ERROR: Unknown systematic '{args.systematic}'")
        print(f"Available systematics: {', '.join(SYSTEMATICS.keys())}")
        return 1
    
    print(f"Running with systematics: {', '.join(systematics_to_run.keys())}")
    
    # Open input file
    fin = ROOT.TFile.Open(args.infile, 'READ')
    if not fin or fin.IsZombie():
        print(f"ERROR: Cannot open input file {args.infile}")
        return 1
    
    # Get trees
    tree_reco = fin.Get('t')
    if not tree_reco:
        print("ERROR: Cannot find tree 't' in input file")
        return 1
    
    tree_gen_before = None
    tree_gen = None
    if args.is_mc:
        tree_gen_before = fin.Get('tgenBefore')
        tree_gen = fin.Get('tgen')
        if not tree_gen_before or not tree_gen:
            print("ERROR: MC mode requested but cannot find 'tgenBefore' or 'tgen' trees")
            return 1
    
    n_events = tree_reco.GetEntries()
    print(f"Processing {n_events} events from {args.infile}")
    print(f"MC mode: {args.is_mc}")
    
    # Create histograms
    print("Creating histograms...")
    
    gen_hists = None
    if args.is_mc:
        gen_hists = create_gen_histograms()
    
    reco_hists = {}
    for syst_name in systematics_to_run.keys():
        reco_hists[syst_name] = create_reco_histograms(syst_name)
    
    # Create response matrices (all systematics in MC)
    response_matrices = {}
    if args.is_mc:
        for syst_name in systematics_to_run.keys():
            response_matrices[syst_name] = create_response_matrices(syst_name)
        print(f"Created response matrices for all {len(systematics_to_run)} systematics")
    
    print(f"Created histograms for {len(systematics_to_run)} systematics")
    
    # Event loop
    print("Starting event loop...")
    n_processed = 0
    
    for iEvt in range(n_events):
        if iEvt % 1000 == 0:
            print(f"  Processing event {iEvt}/{n_events}")
        
        # Process gen level (MC only)
        thrust_gen = None
        if args.is_mc:
            thrust_gen = process_gen_level(
                tree_gen_before, tree_gen, iEvt, gen_hists
            )
            # If gen processing failed (energy cut, etc.), skip this event
            if thrust_gen is None:
                continue
        
        # Process reco level (always)
        process_reco_level(
            tree_reco, iEvt,
            systematics_to_run,
            reco_hists,
            response_matrices,
            thrust_gen=thrust_gen,
            is_mc=args.is_mc
        )
        
        n_processed += 1
    
    print(f"Processed {n_processed} events")
    
    # Write output
    print(f"Writing output to {args.outfile}...")
    fout = ROOT.TFile(args.outfile, 'RECREATE')
    
    # Write gen histograms
    if gen_hists:
        fout.mkdir('gen')
        fout.cd('gen')
        for hist in gen_hists.values():
            hist.Write()
        fout.cd()
    
    # Write reco histograms
    fout.mkdir('reco')
    fout.cd('reco')
    for syst_name, hists in reco_hists.items():
        for hist in hists.values():
            hist.Write()
    fout.cd()
    
    # Write response matrices
    if response_matrices:
        fout.mkdir('response')
        fout.cd('response')
        for syst_name, matrices in response_matrices.items():
            for matrix in matrices.values():
                matrix.Write()
        fout.cd()
    
    fout.Close()
    fin.Close()
    
    print("Done!")
    return 0


if __name__ == "__main__":
    import sys
    sys.exit(main())