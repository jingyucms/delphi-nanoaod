import ROOT
import numpy as np
import sys
import os
import argparse
from array import array

from ROOT import RooUnfoldResponse
from ROOT import RooUnfold
from ROOT import RooUnfoldBayes
from common_functions import *
from binning_and_selections import *

def Proj2D_Y(h,xmin,xmax,hname="XXX"):
    imin=h.GetXaxis().FindBin(xmin)
    imax=h.GetXaxis().FindBin(xmax)-1
    
    proj_y=h.ProjectionY(hname, imin, imax)
    ROOT.SetOwnership(proj_y,True)
    return proj_y

def Proj2D_X(h,ymin,ymax,hname="XXX",Debug=False):
    imin=h.GetYaxis().FindBin(ymin)
    imax=h.GetYaxis().FindBin(ymax)-1

    proj_x=h.ProjectionX(hname, imin, imax)
    ROOT.SetOwnership(proj_x,True)
    return proj_x

def convert_to_roounfold_response(reco_hist, gen_hist, response_hist, name="response"):
    """
    Much simpler approach using direct constructor!
    """
    
    print(f"Creating RooUnfoldResponse '{name}' directly from histograms")
    
    # Just use the direct constructor - that's it!
    response = ROOT.RooUnfoldResponse(
        reco_hist,      # measured TH1
        gen_hist,       # truth TH1
        response_hist,   # response TH2
        name,           # name
        name,           # title
        False           # do_overflow
    )
    
    return response

if __name__ == '__main__':
    
    parser = argparse.ArgumentParser(description='RooUnfold unfolding script')
    parser.add_argument('--iterations', type=int, default=6,
                        help='Number of iterations for Bayesian unfolding')
    parser.add_argument('--version', type=str, default='v40',
                        help='Version tag for input/output files (e.g., v36, v40)')
    
    args = parser.parse_args()
    
    # Version can be changed here or via command line
    VERSION = args.version

    # Define thrust scenarios to unfold
    # Define thrust scenarios to unfold
    thrust_modes = [
        {
            'name': 'thrust2',
            'response_name': 'response_thrust2',
            'reco_name': 'reco_thrust2',
            'gen_name': 'gen_thrust2',
            'before_name': 'Thrust_before2',
            'dataname': 'ThrustMissPNC2'
        },
        {
            'name': 'thrust_log2',
            'response_name': 'response_thrust_log2',
            'reco_name': 'reco_thrust_log2',
            'gen_name': 'gen_thrust_log2',
            'before_name': 'Thrust_before_log2',
            'dataname': 'ThrustMissPNCLog2'
        },
        {
            'name': 'thrust2_Escheme',
            'response_name': 'response_thrust2_Escheme',
            'reco_name': 'reco_thrust2_Escheme',
            'gen_name': 'gen_thrust2_Escheme',
            'before_name': 'Thrust_before2_Escheme',
            'dataname': 'ThrustMissPNC2_Escheme'
        },
        {
            'name': 'thrust_log2_Escheme',
            'response_name': 'response_thrust_log2_Escheme',
            'reco_name': 'reco_thrust_log2_Escheme',
            'gen_name': 'gen_thrust_log2_Escheme',
            'before_name': 'Thrust_before_log2_Escheme',
            'dataname': 'ThrustMissPNCLog2_Escheme'
        },
        # Charged thrust modes
        {
            'name': 'thrustDelphi_c',
            'response_name': 'response_thrustDelphi_c',
            'reco_name': 'reco_thrustDelphi_c',
            'gen_name': 'gen_thrustDelphi_c',
            'before_name': 'ThrustC_beforeDelphi',
            'dataname': 'ThrustCDelphi'
        },
        {
            'name': 'thrustDelphi_c_Escheme',
            'response_name': 'response_thrustDelphi_c_Escheme',
            'reco_name': 'reco_thrustDelphi_c_Escheme',
            'gen_name': 'gen_thrustDelphi_c_Escheme',
            'before_name': 'ThrustC_beforeDelphi_Escheme',
            'dataname': 'ThrustCDelphi_Escheme'
        }
    ]

    # Define all MC generator configurations with their corresponding data files
    # Filename pattern: response_thrust_{generator}_{dataset}_{VERSION}.root
    mc_configs = [
        # 94c data configurations
        {
            'name': 'qqpy_94c',
            'response_file': f'response_thrust_kk2f4146_qqpy_91.25_94c_{VERSION}.root',
            'data_file': f'h_94c_{VERSION}.root',
            'dataset': '94c'
        },
        {
            'name': 'qqardcy_94c',
            'response_file': f'response_thrust_kk2f4146_qqardcy_91.25_94c_{VERSION}.root',
            'data_file': f'h_94c_{VERSION}.root',
            'dataset': '94c'
        },
        {
            'name': 'pythia8_94c',
            'response_file': f'response_thrust_pythia8_94c_{VERSION}.root',
            'data_file': f'h_94c_{VERSION}.root',
            'dataset': '94c'
        },
        {
            'name': 'pythia8_dire_94c',
            'response_file': f'response_thrust_pythia8_dire_94c_{VERSION}.root',
            'data_file': f'h_94c_{VERSION}.root',
            'dataset': '94c'
        },
        # 95d data configurations
        {
            'name': 'qqpy_95d',
            'response_file': f'response_thrust_kk2f4146_qqpy_91.25_95d_{VERSION}.root',
            'data_file': f'h_95d_{VERSION}.root',
            'dataset': '95d'
        },
        {
            'name': 'pythia8_95d',
            'response_file': f'response_thrust_pythia8_95d_{VERSION}.root',
            'data_file': f'h_95d_{VERSION}.root',
            'dataset': '95d'
        },
        {
            'name': 'pythia8_dire_95d',
            'response_file': f'response_thrust_pythia8_dire_95d_{VERSION}.root',
            'data_file': f'h_95d_{VERSION}.root',
            'dataset': '95d'
        }
    ]

    # Single output file for all unfolded results
    filenameout = f"unfolded_data_all_modes_all_generators_{VERSION}_iter{args.iterations}.root"
    fout = ROOT.TFile(filenameout, 'recreate')
    
    print(f"Configuration:")
    print(f"  Version: {VERSION}")
    print(f"  Iterations: {args.iterations}")
    print(f"  Output file: {filenameout}")
    print(f"  Processing {len(thrust_modes)} thrust modes x {len(mc_configs)} MC generators = {len(thrust_modes) * len(mc_configs)} combinations...\n")

    # Track which data files we've already written for each mode
    data_written = {}  # key: (dataset, mode_name)

    # Loop over all thrust modes
    for mode in thrust_modes:
        mode_name = mode['name']
        response_name = mode['response_name']
        reco_name = mode['reco_name']
        gen_name = mode['gen_name']
        before_name = mode['before_name']
        dataname = mode['dataname']
        
        print(f"\n{'#'*70}")
        print(f"# THRUST MODE: {mode_name}")
        print(f"{'#'*70}")
        
        # Loop over all MC configurations
        for config in mc_configs:
            mc_name = config['name']
            filenamein = config['response_file']
            datafile = config['data_file']
            dataset = config['dataset']
            
            # Create unique identifier for this combination
            combo_name = f"{mode_name}_{mc_name}"
            data_key = (dataset, mode_name)
            
            print(f"\n{'='*60}")
            print(f"Processing: {combo_name}")
            print(f"  Response file: {filenamein}")
            print(f"  Data file: {datafile}")
            print(f"  Dataset: {dataset}")
            print(f"{'='*60}")
            
            # Check if files exist
            if not os.path.exists(filenamein):
                print(f"WARNING: Response file {filenamein} not found. Skipping...")
                continue
            
            if not os.path.exists(datafile):
                print(f"WARNING: Data file {datafile} not found. Skipping...")
                continue
            
            # Open response file first to check if this mode exists
            fin = ROOT.TFile.Open(filenamein, 'r')
            
            # Check if the required histograms exist in this response file
            if not fin.Get(response_name):
                print(f"WARNING: {response_name} not found in {filenamein}. Skipping...")
                fin.Close()
                continue
            
            if not fin.Get(reco_name):
                print(f"WARNING: {reco_name} not found in {filenamein}. Skipping...")
                fin.Close()
                continue
                
            if not fin.Get(gen_name):
                print(f"WARNING: {gen_name} not found in {filenamein}. Skipping...")
                fin.Close()
                continue
            
            if not fin.Get(before_name):
                print(f"WARNING: {before_name} not found in {filenamein}. Skipping...")
                fin.Close()
                continue
            
            # Get data
            fdata = ROOT.TFile.Open(datafile, 'r')
            
            # Check if data histogram exists
            if not fdata.Get(dataname):
                print(f"WARNING: {dataname} not found in {datafile}. Skipping...")
                fdata.Close()
                fin.Close()
                continue
            
            data = fdata.Get(dataname).Clone(f"data_{dataset}_{mode_name}")
            try:
                n = fdata.Get('counter').GetBinContent(2)
                data_counter = fdata.Get("counter").Clone(f"NData_{dataset}_{mode_name}")
            except:
                n = fdata.Get('N').GetBinContent(2)
                data_counter = fdata.Get("N").Clone(f"NData_{dataset}_{mode_name}")
            
            # Write data and counter only once per dataset+mode combination
            if data_key not in data_written:
                fout.cd()
                data_counter.Write(f"NData_{dataset}_{mode_name}")
                data.Write(f"data_{dataset}_{mode_name}")
                data_written[data_key] = True
                print(f"  Wrote data for dataset {dataset}, mode {mode_name}")
            
            # Get histograms from response file
            counter = fin.Get("counter").Clone(f"N_{combo_name}")
            _response = fin.Get(response_name)
            reco = fin.Get(reco_name)
            gen = fin.Get(gen_name).Clone(f"gen_{combo_name}")
            before = fin.Get(before_name).Clone(f"before_{combo_name}")
            
            # Get normalization factors from counter histogram
            n_before = counter.GetBinContent(1)  # Total generated events (before selections)
            n_after = counter.GetBinContent(2)   # Events passing selections
            
            print(f"  Events before selection: {n_before}")
            print(f"  Events after selection: {n_after}")
            print(f"  Selection efficiency: {n_after/n_before if n_before > 0 else 0:.4f}")
            
            # Calculate correction factor: corr = (before/n_before) / (gen/n_after)
            # This gives the full phase space correction accounting for selection efficiency
            before_normalized = before.Clone(f"before_norm_temp")
            gen_normalized = gen.Clone(f"gen_norm_temp")
            
            if n_before > 0:
                before_normalized.Scale(1.0 / n_before)
            if n_after > 0:
                gen_normalized.Scale(1.0 / n_after)
            
            corr = before_normalized.Clone(f"corr_{combo_name}")
            corr.Divide(gen_normalized)
            
            # Create response object
            response = convert_to_roounfold_response(reco, gen, _response)
            
            # Check response matrix properties
            RESPONSE = response.Mresponse()
            singular = ROOT.TDecompSVD(RESPONSE)
            print(f"Response matrix singular values for {combo_name}:")
            singular.GetSig().Print()
            
            print(f"Data integral: {data.Integral()}")
            print(f"Before integral: {before.Integral()}")
            print(f"Gen integral: {gen.Integral()}")
            
            # Perform unfolding
            print(f"Performing Bayesian unfolding with {args.iterations} iterations...")
            unfold = ROOT.RooUnfoldBayes(response, data, args.iterations)
            
            hUnf = unfold.Hunfold().Clone(f"unfolded_{combo_name}")
            print(f"Unfolded integral: {hUnf.Integral()}")
            hErr = unfold.Eunfold().Clone(f"unfolding_errors_{combo_name}")
            
            # Write to output file
            fout.cd()
            counter.Write(f"N_{combo_name}")
            _response.Write(f"response_{combo_name}")
            reco.Write(f"reco_{combo_name}")
            gen.Write(f"gen_{combo_name}") 
            before.Write(f"before_{combo_name}")
            corr.Write(f"corr_{combo_name}")
            hUnf.Write(f"unfolded_{combo_name}")
            hErr.Write(f"unfolding_errors_{combo_name}")
            
            # Close input files
            fin.Close()
            fdata.Close()
            
            print(f"Completed unfolding for {combo_name}")

    # Close output file
    fout.Close()

    print(f"\n{'#'*70}")
    print(f"All unfolding completed successfully!")
    print(f"Output saved to: {filenameout}")
    print(f"{'#'*70}")