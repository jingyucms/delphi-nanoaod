#!/usr/bin/env python3
"""
Cross-unfolding script for thrust analysis
-------------------------------------------
Unfold 94c data with 95d response matrix and vice versa.
Only for kk2f4146_qqpy generator, nominal systematic only, all thrust modes.
"""

import ROOT
import sys
import os
import argparse

from ROOT import RooUnfoldResponse
from ROOT import RooUnfoldBayes

def convert_to_roounfold_response(reco_hist, gen_hist, response_hist, name="response"):
    """Convert TH1 and TH2 histograms to RooUnfoldResponse object."""
    # Project response matrix to get distributions
    gen_from_response = response_hist.ProjectionY(f"{name}_gen_projection")
    reco_from_response = response_hist.ProjectionX(f"{name}_reco_projection")
    
    response = ROOT.RooUnfoldResponse(
        reco_from_response,
        gen_from_response,
        response_hist,
        name,
        name,
        False
    )
    return response

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='Cross-unfolding: 94c data with 95d MC and vice versa')
    parser.add_argument('--iterations', type=int, default=6,
                        help='Number of iterations for Bayesian unfolding')
    parser.add_argument('--version', type=str, default='v40',
                        help='Version tag for input files (e.g., v40)')
    
    args = parser.parse_args()
    VERSION = args.version
    SYST = 'nominal'  # Only nominal systematic
    
    # Define thrust modes (all 8 modes)
    thrust_modes = [
        {
            'name': 'thrust2',
            'response_base': 'response_thrust2',
            'reco_base': 'ThrustMissPNC2',
            'gen_name': 'Thrust_before2',
            'data_name': 'ThrustMissPNC2_nominal'
        },
        {
            'name': 'thrust_log2',
            'response_base': 'response_thrust_log2',
            'reco_base': 'ThrustMissPNCLog2',
            'gen_name': 'Thrust_before_log2',
            'data_name': 'ThrustMissPNCLog2_nominal'
        },
        {
            'name': 'thrust2_Escheme',
            'response_base': 'response_thrust2_Escheme',
            'reco_base': 'ThrustMissPNC2_Escheme',
            'gen_name': 'Thrust_before2_Escheme',
            'data_name': 'ThrustMissPNC2_Escheme_nominal'
        },
        {
            'name': 'thrust_log2_Escheme',
            'response_base': 'response_thrust_log2_Escheme',
            'reco_base': 'ThrustMissPNCLog2_Escheme',
            'gen_name': 'Thrust_before_log2_Escheme',
            'data_name': 'ThrustMissPNCLog2_Escheme_nominal'
        },
        {
            'name': 'thrustC2',
            'response_base': 'response_thrustC2',
            'reco_base': 'ThrustC2',
            'gen_name': 'ThrustC_before2',
            'data_name': 'ThrustC2_nominal'
        },
        {
            'name': 'thrustC_log2',
            'response_base': 'response_thrustC_log2',
            'reco_base': 'ThrustCLog2',
            'gen_name': 'ThrustC_before_log2',
            'data_name': 'ThrustCLog2_nominal'
        },
        {
            'name': 'thrustC2_Escheme',
            'response_base': 'response_thrustC2_Escheme',
            'reco_base': 'ThrustC2_Escheme',
            'gen_name': 'ThrustC_before2_Escheme',
            'data_name': 'ThrustC2_Escheme_nominal'
        },
        {
            'name': 'thrustC_log2_Escheme',
            'response_base': 'response_thrustC_log2_Escheme',
            'reco_base': 'ThrustCLog2_Escheme',
            'gen_name': 'ThrustC_before_log2_Escheme',
            'data_name': 'ThrustCLog2_Escheme_nominal'
        }
    ]
    
    # Cross-unfolding configurations
    cross_configs = [
        {
            'name': '94c_data_with_95d_response',
            'data_file': f'thrust_data_94c_{VERSION}.root',
            'mc_file': f'thrust_kk2f4146_qqpy_91.25_95d_{VERSION}.root',
            'data_dataset': '94c',
            'mc_dataset': '95d'
        },
        {
            'name': '95d_data_with_94c_response',
            'data_file': f'thrust_data_95d_{VERSION}.root',
            'mc_file': f'thrust_kk2f4146_qqpy_91.25_94c_{VERSION}.root',
            'data_dataset': '95d',
            'mc_dataset': '94c'
        }
    ]
    
    # Output file
    filenameout = f"unfolded_thrust_cross_{VERSION}_iter{args.iterations}.root"
    fout = ROOT.TFile(filenameout, 'recreate')
    
    print(f"{'='*70}")
    print(f"CROSS-UNFOLDING CONFIGURATION")
    print(f"{'='*70}")
    print(f"  Version: {VERSION}")
    print(f"  Iterations: {args.iterations}")
    print(f"  Generator: kk2f4146_qqpy only")
    print(f"  Systematic: nominal only")
    print(f"  Thrust modes: {len(thrust_modes)}")
    print(f"  Cross configs: {len(cross_configs)}")
    print(f"  Output: {filenameout}")
    print(f"{'='*70}\n")
    
    # Loop over cross-unfolding configurations
    for config in cross_configs:
        config_name = config['name']
        data_file = config['data_file']
        mc_file = config['mc_file']
        data_dataset = config['data_dataset']
        mc_dataset = config['mc_dataset']
        
        print(f"\n{'#'*70}")
        print(f"# CONFIG: {config_name}")
        print(f"{'#'*70}")
        print(f"  Data: {data_file}")
        print(f"  Response: {mc_file}")
        
        # Check files exist
        if not os.path.exists(data_file):
            print(f"ERROR: Data file not found: {data_file}")
            continue
        
        if not os.path.exists(mc_file):
            print(f"ERROR: MC file not found: {mc_file}")
            continue
        
        # Open files
        fdata = ROOT.TFile.Open(data_file, 'READ')
        fmc = ROOT.TFile.Open(mc_file, 'READ')
        
        # Loop over thrust modes
        for mode in thrust_modes:
            mode_name = mode['name']
            response_base = mode['response_base']
            reco_base = mode['reco_base']
            gen_name = mode['gen_name']
            data_name = mode['data_name']
            
            print(f"\n  Mode: {mode_name}")
            
            # Build unique name for this combination
            combo_name = f"{mode_name}_qqpy_{config_name}"
            
            # Get data histogram
            data_hist_path = f"reco/{data_name}"
            data_hist = fdata.Get(data_hist_path)
            
            if not data_hist:
                print(f"    WARNING: Data histogram not found: {data_hist_path}")
                continue
            
            # Get response matrix from MC
            response_name = f"response/{response_base}_{SYST}"
            response_hist = fmc.Get(response_name)
            
            if not response_hist:
                print(f"    WARNING: Response matrix not found: {response_name}")
                continue
            
            # Get reco histogram from MC (for RooUnfoldResponse)
            reco_name = f"reco/{reco_base}_{SYST}"
            reco_hist = fmc.Get(reco_name)
            
            if not reco_hist:
                print(f"    WARNING: Reco histogram not found: {reco_name}")
                continue
            
            # Get gen histogram from MC
            gen_hist_path = f"gen/{gen_name}"
            gen_hist = fmc.Get(gen_hist_path)
            
            if not gen_hist:
                print(f"    WARNING: Gen histogram not found: {gen_hist_path}")
                continue
            
            # Clone histograms
            data = data_hist.Clone(f"data_{combo_name}")
            response_2d = response_hist.Clone(f"response_{combo_name}")
            reco = reco_hist.Clone(f"reco_{combo_name}")
            gen = gen_hist.Clone(f"gen_{combo_name}")
            
            print(f"    Data integral: {data.Integral():.1f}")
            print(f"    Response integral: {response_2d.Integral():.1f}")
            print(f"    Reco (MC) integral: {reco.Integral():.1f}")
            print(f"    Gen (MC) integral: {gen.Integral():.1f}")
            
            # Create RooUnfoldResponse
            response = convert_to_roounfold_response(reco, gen, response_2d, 
                                                      name=f"response_{combo_name}")
            
            # Perform unfolding
            print(f"    Unfolding with {args.iterations} iterations...")
            data_clone = data.Clone(f"data_clone_{combo_name}")
            unfold = ROOT.RooUnfoldBayes(response, data_clone, args.iterations)
            
            hUnf = unfold.Hunfold().Clone(f"unfolded_{combo_name}")
            print(f"    Unfolded integral: {hUnf.Integral():.1f}")
            
            # Get error histogram
            hErr = unfold.Eunfold(ROOT.RooUnfold.kErrors).Clone(f"errors_{combo_name}")
            
            # Write to output
            fout.cd()
            data.Write(f"data_{combo_name}")
            response_2d.Write(f"response_{combo_name}")
            reco.Write(f"reco_{combo_name}")
            gen.Write(f"gen_{combo_name}")
            hUnf.Write(f"unfolded_{combo_name}")
            hErr.Write(f"errors_{combo_name}")
            
            print(f"    ✓ Completed: {combo_name}")
        
        # Close input files
        fdata.Close()
        fmc.Close()
    
    # Close output
    fout.Close()
    
    print(f"\n{'='*70}")
    print(f"Cross-unfolding completed!")
    print(f"Output: {filenameout}")
    print(f"{'='*70}")