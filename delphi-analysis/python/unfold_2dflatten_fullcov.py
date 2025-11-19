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

def Proj2D_X(h,ymin,ymax,hname="XXX",Debug=False):
    # project 2D histogram into 1D along X
    imin=h.GetYaxis().FindBin(ymin)
    imax=h.GetYaxis().FindBin(ymax)-1

    proj_x=h.ProjectionX(hname, imin, imax)
    ROOT.SetOwnership(proj_x,True)

    return proj_x

def flatten_2d_histogram_with_overflow(hist_2d, name_suffix="_flattened"):
    """
    Convert 2D histogram to 1D by flattening, including overflow/underflow.
    """
    
    # Include overflow/underflow bins
    nx = hist_2d.GetNbinsX() + 2  # +2 for under/overflow
    ny = hist_2d.GetNbinsY() + 2  # +2 for under/overflow
    i_start, i_end = 0, nx
    j_start, j_end = 0, ny
    print(f"Including overflow: {nx}x{ny} bins (with under/overflow)")
    
    total_bins = nx * ny
    
    print(f"Flattening 2D histogram -> 1D with {total_bins} bins")
    
    # Create 1D histogram
    flat_name = hist_2d.GetName() + name_suffix
    hist_1d = ROOT.TH1D(flat_name, hist_2d.GetTitle(), total_bins, 0.5, total_bins + 0.5)
    
    # Extract bin edges for overflow bins
    # Get nominal bin edges first
    x_nominal_edges = []
    y_nominal_edges = []
    
    for i in range(hist_2d.GetNbinsX() + 1):
        x_nominal_edges.append(hist_2d.GetXaxis().GetBinLowEdge(i + 1))
    
    for j in range(hist_2d.GetNbinsY() + 1):
        y_nominal_edges.append(hist_2d.GetYaxis().GetBinLowEdge(j + 1))
    
    # Create extended edges: [underflow_edge] + [nominal_edges] + [overflow_edge]
    x_first_width = x_nominal_edges[1] - x_nominal_edges[0]
    x_last_width = x_nominal_edges[-1] - x_nominal_edges[-2]
    y_first_width = y_nominal_edges[1] - y_nominal_edges[0]
    y_last_width = y_nominal_edges[-1] - y_nominal_edges[-2]
    
    # Build extended edge arrays
    x_edges = [x_nominal_edges[0] - x_first_width] + x_nominal_edges + [x_nominal_edges[-1] + x_last_width]
    y_edges = [y_nominal_edges[0] - y_first_width] + y_nominal_edges + [y_nominal_edges[-1] + y_last_width]
    
    print(f"Extended X edges: {len(x_edges)} total ({len(x_nominal_edges)} nominal + 2 overflow)")
    print(f"Extended Y edges: {len(y_edges)} total ({len(y_nominal_edges)} nominal + 2 overflow)")
    
    # Enhanced mapping information
    mapping = {
        'nx': nx,
        'ny': ny,
        'include_overflow': True,
        'x_edges': x_edges,
        'y_edges': y_edges,
        'x_min': hist_2d.GetXaxis().GetXmin(),
        'x_max': hist_2d.GetXaxis().GetXmax(),
        'y_min': hist_2d.GetYaxis().GetXmin(),
        'y_max': hist_2d.GetYaxis().GetXmax()
    }
    
    # Fill 1D histogram
    for i in range(i_start, i_end):
        for j in range(j_start, j_end):
            content = hist_2d.GetBinContent(i, j)
            error = hist_2d.GetBinError(i, j)
            
            # Calculate flat bin index
            flat_bin = i * ny + j + 1  # i and j start from 0
            
            hist_1d.SetBinContent(flat_bin, content)
            hist_1d.SetBinError(flat_bin, error)
    
    print(f"Original integral: {hist_2d.Integral():.6f}")
    print(f"Flattened integral: {hist_1d.Integral():.6f}")
    
    # Check overflow/underflow content
    overflow_content = 0.0
    for i in [0, nx-1]:  # First and last i (underflow/overflow)
        for j in range(ny):
            overflow_content += hist_1d.GetBinContent(i * ny + j + 1)
    for j in [0, ny-1]:  # First and last j (underflow/overflow) 
        for i in range(1, nx-1):  # Avoid double counting corners
            overflow_content += hist_1d.GetBinContent(i * ny + j + 1)
    
    print(f"Overflow/underflow content: {overflow_content:.6f}")
    
    return hist_1d, mapping

def unflatten_to_2d_nominal(hist_1d, mapping, name_suffix="_unflattened"):
    """
    Convert flattened 1D histogram back to 2D, keeping ONLY nominal range bins.
    Overflow/underflow bins are captured in the response matrix but ignored in final result.
    """
    
    nx = mapping['nx']
    ny = mapping['ny']
    
    # Original histogram dimensions (without overflow)
    nx_nominal = nx - 2
    ny_nominal = ny - 2
    print(f"Unflattening with overflow -> nominal 2D {nx_nominal}x{ny_nominal} (ignoring overflow bins)")
    
    # Create 2D histogram with ONLY nominal range
    unflatten_name = hist_1d.GetName() + name_suffix
    
    x_edges_array = array('d', mapping['x_edges'])
    y_edges_array = array('d', mapping['y_edges'])
    
    if len(x_edges_array) != nx + 1 or len(y_edges_array) != ny + 1:
        print(f"WARNING: Edge array size mismatch. X: {len(x_edges_array)} vs {nx+1}, Y: {len(y_edges_array)} vs {ny+1}")
    
    # Extract nominal edges: skip first (underflow) and last (overflow)
    x_nominal_edges = array('d', x_edges_array[1:nx_nominal+2])  # [1:nx-1] = nominal edges
    y_nominal_edges = array('d', y_edges_array[1:ny_nominal+2])  # [1:ny-1] = nominal edges
    
    print(f"Extracted {len(x_nominal_edges)-1} nominal X bins from {len(x_edges_array)-1} total")
    print(f"Extracted {len(y_nominal_edges)-1} nominal Y bins from {len(y_edges_array)-1} total")
    print(f"X range: [{x_nominal_edges[0]:.6f}, {x_nominal_edges[-1]:.6f}]")
    print(f"Y range: [{y_nominal_edges[0]:.6f}, {y_nominal_edges[-1]:.6f}]")
    
    hist_2d = ROOT.TH2D(
        unflatten_name, hist_1d.GetTitle(),
        nx_nominal, x_nominal_edges,
        ny_nominal, y_nominal_edges
    )
    
    # Fill 2D histogram - ONLY nominal bins
    nominal_content = 0.0
    overflow_content = 0.0
    
    for flat_bin in range(1, hist_1d.GetNbinsX() + 1):
        content = hist_1d.GetBinContent(flat_bin)
        error = hist_1d.GetBinError(flat_bin)
        
        # Reverse the flattening
        linear_index = flat_bin - 1
        
        i = linear_index // ny  # 0-based, includes overflow
        j = linear_index % ny   # 0-based, includes overflow
        
        # Only fill nominal bins (exclude overflow: i=0, i=nx-1, j=0, j=ny-1)
        if i > 0 and i < nx-1 and j > 0 and j < ny-1:
            # Convert to 1-based nominal histogram coordinates
            i_nominal = i         # i=1..nx-2 maps to 1..nx_nominal  
            j_nominal = j         # j=1..ny-2 maps to 1..ny_nominal
            hist_2d.SetBinContent(i_nominal, j_nominal, content)
            hist_2d.SetBinError(i_nominal, j_nominal, error)
            nominal_content += content
        else:
            # This is overflow/underflow content - just track it
            overflow_content += content
    
    print(f"Nominal range integral: {nominal_content:.6f}")
    if overflow_content > 0:
        print(f"Overflow/underflow content ignored: {overflow_content:.6f}")
        print(f"Fraction kept in nominal range: {nominal_content/(nominal_content + overflow_content):.3f}")
    
    print(f"Final 2D histogram integral: {hist_2d.Integral():.6f}")
    
    return hist_2d

def unflatten_errors_to_2d_nominal(error_1d, mapping, name_suffix="_error_2d"):
    """
    Convert flattened 1D error histogram to 2D nominal range.
    """
    
    nx_full = mapping['nx']
    ny_full = mapping['ny']
    nx_nominal = nx_full - 2
    ny_nominal = ny_full - 2
    
    print(f"Unflattening errors: {error_1d.GetNbinsX()} -> {nx_nominal}x{ny_nominal}")
    
    # Create 2D error histogram for nominal range
    error_name = error_1d.GetName() + name_suffix
    x_edges_array = array('d', mapping['x_edges'][1:-1])  # Skip overflow edges
    y_edges_array = array('d', mapping['y_edges'][1:-1])  # Skip overflow edges
    
    error_2d = ROOT.TH2D(
        error_name, "Per-bin Errors (2D)",
        nx_nominal, x_edges_array,
        ny_nominal, y_edges_array
    )
    
    # Map from 1D nominal error bins to 2D nominal bins
    for i in range(1, nx_nominal + 1):
        for j in range(1, ny_nominal + 1):
            # Calculate flat bin in nominal space
            flat_bin = (i - 1) * ny_nominal + j
            
            if flat_bin <= error_1d.GetNbinsX():
                error_value = error_1d.GetBinContent(flat_bin)
                error_2d.SetBinContent(i, j, error_value)
    
    print(f"Error 2D integral after unflattening: {error_2d.Integral():.6f}")
    
    return error_2d

def create_2d_response_matrix_from_thnf_with_overflow(thnf_4d, name_suffix="_response_2d"):
    """
    Convert 4D THnF to 2D TH2 response matrix, including overflow/underflow.
    """
    
    print(f"Converting 4D THnF to 2D response matrix with overflow capture...")
    
    if thnf_4d.GetNdimensions() != 4:
        raise ValueError(f"Expected 4D THnF, got {thnf_4d.GetNdimensions()}D")
    
    # Get dimensions
    axis_reco_x = thnf_4d.GetAxis(0)
    axis_reco_y = thnf_4d.GetAxis(1)
    axis_gen_x = thnf_4d.GetAxis(2)
    axis_gen_y = thnf_4d.GetAxis(3)
    
    # Include overflow/underflow bins
    nx_reco = axis_reco_x.GetNbins() + 2
    ny_reco = axis_reco_y.GetNbins() + 2
    nx_gen = axis_gen_x.GetNbins() + 2
    ny_gen = axis_gen_y.GetNbins() + 2
    
    reco_start, reco_end = 0, nx_reco
    gen_start, gen_end = 0, nx_gen
    reco_j_start, reco_j_end = 0, ny_reco
    gen_j_start, gen_j_end = 0, ny_gen
    
    total_reco_bins = nx_reco * ny_reco
    total_gen_bins = nx_gen * ny_gen
    
    print(f"THnF dimensions: reco({nx_reco}x{ny_reco}) -> gen({nx_gen}x{ny_gen})")
    print(f"Flattened response matrix: {total_reco_bins} x {total_gen_bins}")
    
    # Create 2D response matrix
    response_name = thnf_4d.GetName() + name_suffix
    response_2d = ROOT.TH2D(
        response_name, "2D Response Matrix",
        total_reco_bins, 0.5, total_reco_bins + 0.5,  # Reco (X-axis)
        total_gen_bins, 0.5, total_gen_bins + 0.5      # Gen (Y-axis)
    )
    
    print("Filling 2D response matrix...")
    
    filled_bins = 0
    coords = array('i', [0, 0, 0, 0])
    total_content = 0.0
    overflow_content = 0.0
    
    # Iterate through all bins (including overflow)
    for i_reco in range(reco_start, reco_end):
        for j_reco in range(reco_j_start, reco_j_end):
            # Calculate flat reco bin
            flat_reco_bin = i_reco * ny_reco + j_reco + 1
            
            for i_gen in range(gen_start, gen_end):
                for j_gen in range(gen_j_start, gen_j_end):
                    # Calculate flat gen bin
                    flat_gen_bin = i_gen * ny_gen + j_gen + 1
                    
                    # Map 0-based indices to THnF overflow convention
                    if i_reco == 0:
                        coords[0] = 0  # Underflow
                    elif i_reco == nx_reco - 1:
                        coords[0] = axis_reco_x.GetNbins() + 1  # Overflow
                    else:
                        coords[0] = i_reco  # Normal bin
                    
                    if j_reco == 0:
                        coords[1] = 0
                    elif j_reco == ny_reco - 1:
                        coords[1] = axis_reco_y.GetNbins() + 1
                    else:
                        coords[1] = j_reco
                        
                    if i_gen == 0:
                        coords[2] = 0
                    elif i_gen == nx_gen - 1:
                        coords[2] = axis_gen_x.GetNbins() + 1
                    else:
                        coords[2] = i_gen
                        
                    if j_gen == 0:
                        coords[3] = 0
                    elif j_gen == ny_gen - 1:
                        coords[3] = axis_gen_y.GetNbins() + 1
                    else:
                        coords[3] = j_gen
                    
                    try:
                        content = thnf_4d.GetBinContent(coords)
                        if content > 0:
                            response_2d.SetBinContent(flat_reco_bin, flat_gen_bin, content)
                            filled_bins += 1
                            total_content += content
                            
                            # Track overflow content
                            is_overflow = (i_reco == 0 or i_reco == nx_reco-1 or 
                                         j_reco == 0 or j_reco == ny_reco-1 or
                                         i_gen == 0 or i_gen == nx_gen-1 or
                                         j_gen == 0 or j_gen == ny_gen-1)
                            if is_overflow:
                                overflow_content += content
                        else:
                            response_2d.SetBinContent(flat_reco_bin, flat_gen_bin, 0)
                    except:
                        continue
        
        if i_reco % 20 == 0:
            print(f"  Processed reco slice {i_reco}/{reco_end - 1}")
    
    print(f"Filled {filled_bins} response matrix elements")
    print(f"Total response content: {total_content:.6f}")
    print(f"Overflow/underflow content: {overflow_content:.6f} ({100*overflow_content/total_content:.1f}%)")
    
    # Create mapping information
    mapping = {
        'reco_nx': nx_reco, 'reco_ny': ny_reco,
        'gen_nx': nx_gen, 'gen_ny': ny_gen,
        'include_overflow': True,
        'reco_x_min': axis_reco_x.GetBinLowEdge(1),
        'reco_x_max': axis_reco_x.GetBinUpEdge(axis_reco_x.GetNbins()),
        'reco_y_min': axis_reco_y.GetBinLowEdge(1),
        'reco_y_max': axis_reco_y.GetBinUpEdge(axis_reco_y.GetNbins()),
        'gen_x_min': axis_gen_x.GetBinLowEdge(1),
        'gen_x_max': axis_gen_x.GetBinUpEdge(axis_gen_x.GetNbins()),
        'gen_y_min': axis_gen_y.GetBinLowEdge(1),
        'gen_y_max': axis_gen_y.GetBinUpEdge(axis_gen_y.GetNbins())
    }
    
    return response_2d, mapping

def extract_covariance_and_errors(unfolder, gen_mapping):
    """
    Extract covariance matrix and errors, handling overflow and extracting nominal range.
    Returns nominal-range covariance and errors.
    """
    
    try:
        cov_obj = unfolder.Eunfold()
        
        if cov_obj is None:
            print("No covariance matrix returned")
            return None, None, 'none'
        
        # Convert covariance object to 2D histogram
        cov_hist_full = None
        
        # If it's already a 2D histogram, use as-is
        if hasattr(cov_obj, 'GetNbinsX') and hasattr(cov_obj, 'GetNbinsY'):
            print(f"Covariance is 2D histogram: {cov_obj.GetNbinsX()}x{cov_obj.GetNbinsY()}")
            cov_hist_full = cov_obj
        
        # If it's a matrix, convert to 2D histogram
        elif hasattr(cov_obj, 'GetNrows'):
            print(f"Converting matrix to 2D histogram: {cov_obj.GetNrows()}x{cov_obj.GetNcols()}")
            
            nrows = cov_obj.GetNrows()
            ncols = cov_obj.GetNcols()
            
            cov_hist_full = ROOT.TH2D(
                "covariance_matrix_full", "Full Covariance Matrix",
                nrows, 0.5, nrows + 0.5,
                ncols, 0.5, ncols + 0.5
            )
            
            # Copy matrix values to histogram
            for i in range(nrows):
                for j in range(ncols):
                    cov_value = cov_obj(i, j)
                    cov_hist_full.SetBinContent(i + 1, j + 1, cov_value)
            
            print(f"Created full covariance histogram: {cov_hist_full.GetNbinsX()}x{cov_hist_full.GetNbinsY()}")
        
        # If it's a 1D histogram (diagonal only), convert to 2D diagonal
        elif hasattr(cov_obj, 'GetNbinsX'):
            print(f"Converting 1D diagonal to 2D: {cov_obj.GetNbinsX()} bins")
            
            nbins = cov_obj.GetNbinsX()
            
            cov_hist_full = ROOT.TH2D(
                "covariance_diagonal_full", "Full Diagonal Covariance Matrix",
                nbins, 0.5, nbins + 0.5,
                nbins, 0.5, nbins + 0.5
            )
            
            # Fill diagonal elements only
            for i in range(1, nbins + 1):
                error = cov_obj.GetBinContent(i)
                variance = error * error  # Convert error to variance
                cov_hist_full.SetBinContent(i, i, variance)
            
            print(f"Created full diagonal covariance: {cov_hist_full.GetNbinsX()}x{cov_hist_full.GetNbinsY()}")
        
        else:
            print(f"Unknown covariance type: {type(cov_obj)}")
            return None, None, 'unknown'
        
        if cov_hist_full is None:
            return None, None, 'failed'
        
        # Extract diagonal errors from full matrix and convert to nominal range
        print("Extracting errors and nominal covariance...")
        
        # Get dimensions
        nx_full = gen_mapping['nx']  # Includes overflow
        ny_full = gen_mapping['ny']  # Includes overflow
        nx_nominal = nx_full - 2     # Nominal bins only
        ny_nominal = ny_full - 2     # Nominal bins only
        
        total_nominal_bins = nx_nominal * ny_nominal
        
        print(f"Converting from full ({nx_full}x{ny_full}) to nominal ({nx_nominal}x{ny_nominal})")
        
        # Create nominal covariance matrix
        cov_hist_nominal = ROOT.TH2D(
            "covariance_nominal", "Nominal Covariance Matrix",
            total_nominal_bins, 0.5, total_nominal_bins + 0.5,
            total_nominal_bins, 0.5, total_nominal_bins + 0.5
        )
        
        # Create nominal error histogram
        error_hist_nominal = ROOT.TH1D(
            "errors_nominal", "Nominal Per-bin Errors",
            total_nominal_bins, 0.5, total_nominal_bins + 0.5
        )
        
        # Map nominal bins from full matrix to nominal matrix
        nominal_bin_map = []  # Maps nominal_flat_bin -> full_flat_bin
        
        for i in range(1, nx_nominal + 1):  # 1-based nominal bins
            for j in range(1, ny_nominal + 1):
                # In the full matrix, nominal bins start at i=1, j=1 (after underflow i=0, j=0)
                i_full = i + 1  # Skip underflow bin
                j_full = j + 1  # Skip underflow bin
                
                # Convert to flat bin indices
                full_flat_bin = (i_full - 1) * ny_full + j_full  # 0-based then +1 for ROOT
                nominal_bin_map.append(full_flat_bin)
        
        print(f"Created mapping for {len(nominal_bin_map)} nominal bins")
        
        # Extract the nominal submatrix and errors
        extracted_elements = 0
        total_variance = 0.0
        
        for nom_i, full_i in enumerate(nominal_bin_map):
            for nom_j, full_j in enumerate(nominal_bin_map):
                # Get covariance element from full matrix
                cov_value = cov_hist_full.GetBinContent(full_i + 1, full_j + 1)  # ROOT is 1-based
                
                # Set in nominal matrix
                cov_hist_nominal.SetBinContent(nom_i + 1, nom_j + 1, cov_value)
                
                if nom_i == nom_j:  # Diagonal element
                    total_variance += cov_value
                    # Extract error (square root of variance)
                    error = np.sqrt(np.maximum(0.0, cov_value))
                    error_hist_nominal.SetBinContent(nom_i + 1, error)
                
                if cov_value != 0:
                    extracted_elements += 1
        
        print(f"Extracted {extracted_elements} non-zero covariance elements")
        print(f"Total variance (diagonal sum): {total_variance:.6f}")
        print(f"Error histogram integral: {error_hist_nominal.Integral():.6f}")
        
        return cov_hist_nominal, error_hist_nominal, 'success'
        
    except Exception as e:
        print(f"Error extracting covariance and errors: {e}")
        return None, None, 'failed'

def perform_2d_unfolding_via_flattening_with_overflow(reco_2d, gen_2d, response_thnf, data_2d, 
                                                      n_iterations=4, method="Bayes", h_covariance=None):
    """
    2D EEC unfolding that includes overflow in response matrix but returns only nominal range.
    """
    
    print("="*60)
    print(f"PERFORMING 2D UNFOLDING WITH OVERFLOW CAPTURE")
    print(f"Final result: nominal range only")
    print("="*60)
    
    # Step 1: Flatten with overflow for response matrix accuracy
    print("\nStep 1: Flattening 2D histograms...")
    reco_1d, reco_mapping = flatten_2d_histogram_with_overflow(reco_2d, "_reco_flat")
    gen_1d, gen_mapping = flatten_2d_histogram_with_overflow(gen_2d, "_gen_flat")
    data_1d, data_mapping = flatten_2d_histogram_with_overflow(data_2d, "_data_flat")
    
    # Step 2: Create response matrix with overflow to capture all migrations
    print(f"\nStep 2: Creating response matrix with overflow capture...")
    response_2d, response_mapping = create_2d_response_matrix_from_thnf_with_overflow(response_thnf)
    
    # Step 3: RooUnfold setup
    print("\nStep 3: Creating RooUnfoldResponse...")
    roounfold_response = ROOT.RooUnfoldResponse(
        reco_1d, gen_1d, response_2d,
        "response_2d_overflow", "2D Response with Overflow Capture", False
    )

    # Step 4: Perform unfolding with full overflow-aware response matrix
    if method == "Bayes":
        print(f"\nStep 4: Performing {method} unfolding with {n_iterations} iterations...")
    
        # Create the unfolder object
        unfolder = ROOT.RooUnfoldBayes(roounfold_response, data_1d, n_iterations)
    
        # Check for and set the covariance matrix
        if h_covariance:
            print(f"✅ Setting combined data+fake covariance matrix: {h_covariance.GetNrows()}x{h_covariance.GetNcols()}")
            unfolder.SetMeasuredCov(total_covariance_tmatrix)
        else:
            print("ℹ️ No data covariance matrix provided. Using Poisson errors from data histogram.")

    elif method == "BinByBin":
        print(f"\nStep 4: Performing {method} unfolding...")
        unfolder = ROOT.RooUnfoldBinByBin(roounfold_response, data_1d)
    elif method == "Invert":
        print(f"\nStep 4: Performing {method} unfolding...")
        unfolder = ROOT.RooUnfoldInvert(roounfold_response, data_1d)
    else:
        raise ValueError(f"Unknown method: {method}")
        
    unfolded_1d = unfolder.Hunfold()
    
    print(f"1D unfolding completed")
    print(f"  Data integral: {data_1d.Integral():.6f}")
    print(f"  Unfolded integral: {unfolded_1d.Integral():.6f}")
    
    # Step 5: Extract covariance and errors in one step
    print("\nStep 5: Extracting covariance matrix and errors...")
    cov_2d, error_1d, cov_status = extract_covariance_and_errors(unfolder, gen_mapping)

    cov_2d = unfolder.Eunfold()
    
    # Convert 1D errors to 2D if available
    error_2d = None
    if error_1d is not None:
        error_2d = unflatten_errors_to_2d_nominal(error_1d, gen_mapping, "_error_2d")
    else:
        print("❌ No errors extracted")
    
    # Step 6: Unflatten to nominal range only (ignore overflow bins)
    print("\nStep 6: Unflattening to nominal range (discarding overflow)...")
    unfolded_2d = unflatten_to_2d_nominal(unfolded_1d, gen_mapping, "_unfolded_2d")
    
    print(f"\n=== UNFOLDING SUMMARY ===")
    print(f"Input data: {data_2d.Integral():.6f}")
    print(f"Unfolded (nominal): {unfolded_2d.Integral():.6f}")
    print(f"Efficiency (kept in nominal): {unfolded_2d.Integral()/data_2d.Integral():.3f}")
    
    return unfolded_2d, error_2d, cov_2d, response_2d

def apply_efficiency_correction_hist(hUnf2d, eff_corr):
    """
    Apply efficiency correction while preserving the relative error pattern from unfolding
    """
    
    print("Applying efficiency correction while preserving relative errors...")
    
    for i in range(1, hUnf2d.GetNbinsX() + 1):
        for j in range(1, hUnf2d.GetNbinsY() + 1):
            # Get current values
            content = hUnf2d.GetBinContent(i, j)
            error = hUnf2d.GetBinError(i, j)
            eff_factor = eff_corr.GetBinContent(i, j)
            
            # Calculate relative error BEFORE correction
            rel_error = error / content if content > 0 else 0
            
            # Apply efficiency correction to content only
            new_content = content * eff_factor
            
            # Calculate new error to preserve the same relative error
            new_error = new_content * rel_error
            
            # Set the corrected values
            hUnf2d.SetBinContent(i, j, new_content)
            hUnf2d.SetBinError(i, j, new_error)
    
    print("✅ Efficiency correction applied with preserved relative errors")

def unfold_2d_eec_with_overflow(reco_2d_hist, gen_2d_hist, response_4d_thnf, data_2d_hist, 
                                n_iterations=4, unfolding_method="Invert", h_covariance=None):
    """
    2D EEC unfolding that includes overflow in response matrix but returns only nominal range.
    
    Returns:
    --------
    unfolded_2d : ROOT.TH2D
        Unfolded result in original 2D coordinates (nominal range only)
    error_2d : ROOT.TH2D  
        Per-bin errors in original 2D coordinates (optional, can be None)
    cov_2d : ROOT.TH2D
        Covariance matrix in flattened coordinate space
    """
    
    try:
        unfolded_2d, error_2d, cov_2d, resp_2d = perform_2d_unfolding_via_flattening_with_overflow(
            reco_2d_hist, gen_2d_hist, response_4d_thnf, data_2d_hist, 
            n_iterations, unfolding_method, h_covariance
        )
    
        return unfolded_2d, error_2d, cov_2d, resp_2d
        
    except Exception as e:
        print(f"ERROR in unfold_2d_eec_with_overflow: {e}")
        return None, None, None, None

def th1d_to_numpy_simple(hist):
    """Convert ROOT TH1D to numpy array"""
    n_bins = hist.GetNbinsX()
    array = np.zeros(n_bins)
    for i in range(1, n_bins + 1):
        array[i-1] = hist.GetBinContent(i)
    return array

def th2d_to_numpy_simple(hist):
    """Convert ROOT TH2D to numpy array"""
    nx = hist.GetNbinsX()
    ny = hist.GetNbinsY()
    array = np.zeros((nx, ny))
    for i in range(1, nx + 1):
        for j in range(1, ny + 1):
            array[i-1, j-1] = hist.GetBinContent(i, j)
    return array

def numpy_to_tmatrixd(numpy_array):
    """
    Converts a 2D NumPy array into a ROOT TMatrixD.

    Args:
        numpy_array (np.ndarray): The 2D NumPy array to convert.

    Returns:
        ROOT.TMatrixD: The converted TMatrixD object.
    """
    if numpy_array.ndim != 2:
        raise ValueError("Input must be a 2D NumPy array.")
    
    rows, cols = numpy_array.shape
    tmatrix = ROOT.TMatrixD(rows, cols)
    
    # Loop over elements and fill the TMatrixD
    # NumPy and TMatrixD are both 0-indexed, so the mapping is direct.
    for i in range(rows):
        for j in range(cols):
            tmatrix[i][j] = numpy_array[i, j]
            
    return tmatrix

if __name__ == '__main__':

    parser = argparse.ArgumentParser(description='2D Unfolding with Overflow Handling')
    parser.add_argument('--method', type=str, default="Bayes",
                       help='Unfolding method to use (default: Bayes)')
    parser.add_argument('--niter', type=int, default=4,
                       help='number of iterations')
    parser.add_argument('--jacobian', type=str, default="r",
                       help='Unfold r or z distribution')
    args = parser.parse_args()

    if args.method == "Bayes":
        surfix = f'_niter{args.niter}'
    else:
        surfix = ''
    
    surfix += '_overflow'  # Always use overflow
    surfix += f'_{args.jacobian}'

    # response v16 -> data v5 with covariance (nominal)
    # response v15 -> data v3 vvv tight angular
    # response v14 -> data v1 angular (nominal)
    # response v13 -> data v2 correspondence table
    #filenamein = 'response_kk2f4146_qqpy_91.25_v16_zcov.root'
    #filenamein = 'response_kk2f4146_qqardcy_91.25_v1.root'
    #filenamein = 'response_pythia8_91.25_v3.root'
    #filenamein = 'response_pythia8_dire_91.25_v1.root'

    filenamein = 'response_ALEPHMC_LEP1MC1994_recons_aftercut.v1.root'

    #filenamein = 'response_kk2f4146_qqpy_91.25_95d_v1.root'
    #filenamein = 'response_pythia8_91.25_95d_v1.root'
    #filenamein = 'response_pythia8_dire_91.25_95d_v1.root'
    
    #datafile = 'h_94c_v8_cov.root'
    #datafile = 'h_95d_v1_cov.root'
    datafile = 'h_aleph_v1.root'
    
    #filenameout = f'unfolded_data_kk2f4146_qqpy_91.25_v5{surfix}.root'
    #filenameout = f'unfolded_data_kk2f4146_qqardcy_91.25_v1{surfix}.root'
    #filenameout = f'unfolded_data_pythia8_91.25_v1{surfix}.root'
    #filenameout = f'unfolded_data_pythia8_dire_91.25_v1{surfix}.root'

    #filenameout = f'unfolded_data_95d_kk2f4146_qqpy_91.25_v1{surfix}.root'
    #filenameout = f'unfolded_data_95d_pythia8_91.25_v1{surfix}.root'
    #filenameout = f'unfolded_data_95d_pythia8_dire_91.25_v1{surfix}.root'

    filenameout = f'unfolded_data_aleph_v1{surfix}.root'
    
    data2dname = f'EEC2d_{args.jacobian}'
    
    fin = ROOT.TFile.Open(filenamein,'r')
    h_response = fin.Get(f"response2d_eij_{args.jacobian}")
    h_reco_matched = fin.Get(f'reco2d_eij_{args.jacobian}_match')
    h_gen_matched = fin.Get(f"gen2d_eij_{args.jacobian}_match")

    h_reco_total_mc = fin.Get(f'reco2d_eij_{args.jacobian}')
    h_gen_total_mc = fin.Get(f'gen2d_eij_{args.jacobian}')

    mc_event_count = fin.Get("counter").GetBinContent(2)


    # Open data file
    fdata = ROOT.TFile.Open(datafile, 'r')
    h_data_raw = fdata.Get(data2dname).Clone('h_data_raw')
    try:
        data_event_count = fdata.Get('counter').GetBinContent(2)
        counter = fdata.Get("counter").Clone("N")
    except:
        data_event_count = fdata.Get('N').GetBinContent(2)
        counter = fdata.Get("N").Clone("N")

    print(f"MC Events: {mc_event_count}, Data Events: {data_event_count}")

    scale_factor = float(data_event_count) / mc_event_count
    print(f"MC-to-Data Scale Factor: {scale_factor}")

    # The absolute number of fakes in the MC is (total reco) - (reco from matched gen)
    h_fakes_mc_abs = h_reco_total_mc.Clone("h_fakes_mc_abs")
    h_fakes_mc_abs.Add(h_reco_matched, -1)

    # Now, scale this absolute MC prediction to the data luminosity
    h_fakes_scaled = h_fakes_mc_abs.Clone("h_fakes_scaled")
    h_fakes_scaled.Scale(scale_factor)

    h_data_subtracted = h_data_raw.Clone("h_data_subtracted")

    # Subtract the scaled, absolute fake prediction
    h_data_subtracted.Add(h_fakes_scaled, -1)

    try: 
        sum_of_eec_products = th2d_to_numpy_simple(fdata.Get(f"covariance_matrix_{args.jacobian}"))
        sum_of_eecs = th1d_to_numpy_simple(fdata.Get(f"mean_eec_{args.jacobian}"))
        numerator = data_event_count * sum_of_eec_products - np.outer(sum_of_eecs, sum_of_eecs)
        data_covariance_matrix = numerator / (data_event_count - 1)

        sum_of_eec_products = th2d_to_numpy_simple(fin.Get(f"covariance_matrix_{args.jacobian}"))
        sum_of_eecs = th1d_to_numpy_simple(fin.Get(f"mean_eec_{args.jacobian}"))
        numerator = mc_event_count * sum_of_eec_products -  np.outer(sum_of_eecs, sum_of_eecs)
        fake_covariance_matrix = numerator / (mc_event_count - 1)
        
        total_covariance_numpy = data_covariance_matrix + fake_covariance_matrix*(scale_factor**2)

        total_covariance_tmatrix = numpy_to_tmatrixd(total_covariance_numpy)
        
    except:
        total_covariance_tmatrix = None
        print("❌ Could not find input data covariance matrix!!!")
    

    print("Performing unfolding...")
    hUnf2d, hErr2d, hCov2d, hResp2d = unfold_2d_eec_with_overflow(
        h_reco_matched, h_gen_matched, h_response, h_data_subtracted,
        n_iterations=args.niter,
        unfolding_method=args.method,
        h_covariance=total_covariance_tmatrix
    )

    h_eff = h_gen_total_mc.Clone("h_eff")
    h_eff.Divide(h_gen_matched)

    # Save output
    fout = ROOT.TFile(filenameout, 'recreate')
    fout.cd()

    hResp2d.Write("response_2d")
    hUnf2d.Write("unfolded_2d_noeff")

    #total_covariance_tmatrix.Write("input_cov")

    apply_efficiency_correction_hist(hUnf2d, h_eff)
    
    eec_reco = Proj2D_X(h_data_raw, eijbins[1], eijbins[-1], f"RECO_EEC")
    eec_gen = Proj2D_X(h_gen_total_mc, eijbins[1], eijbins[-1], f"GEN_EEC")
    eec_unfold = Proj2D_X(hUnf2d, eijbins[1], eijbins[-1], f"UNFOLD_EEC")
    eec_reco.SetDirectory(0)
    eec_gen.SetDirectory(0)
    eec_unfold.SetDirectory(0)

    hUnf2d.Write("unfolded_2d")

    if hErr2d is not None:
        hErr2d.Write("errors_2d")      # Per-bin errors in original 2D space
    
    if hCov2d is not None:
        hCov2d.Write("covariance_2d")  # Covariance matrix in flattened space

    counter.Write("N")
        
    h_gen_total_mc.Write("gen_2d")
    h_reco_total_mc.Write("reco_2d")
    h_data_raw.Write("data_2d")

    h_fakes_scaled.Write("fake_corr")
    h_eff.Write("eff_corr")

    h_gen_matched.Write("gen_2d_match")
    h_reco_matched.Write("reco_2d_match")

    eec_reco.Write()
    eec_gen.Write()
    eec_unfold.Write()

    fout.Close()
    
    print(f"\n=== UNFOLDING COMPLETE ===")
    print(f"Method: {args.method}")
    if args.method == "Bayes":
        print(f"Iterations: {args.niter}")
    print(f"Overflow handling: Always enabled")
    print(f"Output file: {filenameout}")
    print(f"Final unfolded integral: {hUnf2d.Integral():.6f}")
