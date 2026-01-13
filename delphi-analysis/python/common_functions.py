import ROOT
import math
import numpy as np
from array import array
from itertools import combinations
from scipy.optimize import linear_sum_assignment
from scipy.linalg import eigh as scipy_eigh
from dataclasses import dataclass

# --------------------------------------------------------------------
#  dataclass + factory (only **one** public entry point → P4Block.build)
# --------------------------------------------------------------------
@dataclass
class P4Block:
    p       : np.ndarray        # (N,3)  px,py,pz
    q       : np.ndarray        # (N,)   charge
    pt      : np.ndarray        # (N,)
    theta   : np.ndarray        # (N,)   <<—  polar angle
    phi     : np.ndarray        # (N,)
    e       : np.ndarray        # (N,)   energy
    m       : np.ndarray        # (N,)   mass
    idx     : np.ndarray = None # (N,)   original index
    cspidx  : np.ndarray = None # (N,)   correspondence index

    @classmethod
    def build(cls, px, py, pz, q, pt, theta, phi, e, m=None, idx=None, cspidx=None):
        """
        Build P4Block with optional mass, index, and correspondence index.
        
        Parameters:
        -----------
        px, py, pz, q, pt, theta, phi, e : array-like
            Required particle properties
        m : array-like, optional
            Mass array. If None, calculated from 4-momentum
        idx : array-like, optional
            Original particle indices
        cspidx : array-like, optional
            Correspondence indices for matching
        """
        # Calculate mass if not provided
        if m is None:
            m = np.sqrt(np.maximum(0, e**2 - px**2 - py**2 - pz**2))
        
        return cls(np.stack([px,py,pz], axis=1), q, pt, theta, phi, e, m, idx, cspidx)

    @property
    def norm(self): return np.linalg.norm(self.p, axis=1)

def calcBinEdge(low, high, nbins):   
    xlow = np.log10(low)
    xhigh = np.log10(high)
    width = (xhigh-xlow)/nbins

    bins=[]
    for i in range(nbins+1):
        val = pow(10, xlow + i * width)
        bins += [val]
    
    newbins = [2*high - b for b in bins]
    newbins = newbins[::-1]
    del newbins[0]
    bin_edge = bins + newbins

    return bin_edge

def normalizeByBinWidth(h):
    for b in range(h.GetNbinsX()):
        h.SetBinContent(b+1, h.GetBinContent(b+1)/h.GetBinWidth(b+1))
        h.SetBinError(b+1, h.GetBinError(b+1)/h.GetBinWidth(b+1))

    return h

def calcAngle(n1, n2):
    cos_theta = np.dot(n1, n2) / (np.linalg.norm(n1) * np.linalg.norm(n2))
    cos_theta = np.clip(cos_theta, -1.0, 1.0)  
    
    theta = np.arccos(cos_theta)
    return theta

def missing_p(p3: np.ndarray) -> np.ndarray:
    """
    Missing 3-momentum vector  (−Σ p_i).

    Parameters
    ----------
    p3 : (n, 3) ndarray
        Each row = (px, py, pz) of one particle (any units).

    Returns
    -------
    (3,) ndarray
        [-Σpx,  -Σpy,  -Σpz]   ─ i.e. the negative of the vector sum
        of all particle momenta.
    """
    # Sum over all rows, keep all three components
    total_p = p3.sum(axis=0)          # shape (3,)
    return -total_p                   # missing-momentum (px,py,pz)

def polar_angle(p3: np.ndarray) -> float:
    """
    Calculate the polar angle (theta) of a 3-momentum vector.

    Parameters
    ----------
    p3 : (3,) ndarray
        (px, py, pz) of a single momentum vector.

    Returns
    -------
    float
        Polar angle θ (in radians) defined as arccos(pz / |p|),
        where |p| = sqrt(px² + py² + pz²).
        Range: [0, π]
    """
    # Calculate polar angle: θ = arccos(pz / |p|)
    p_magnitude = np.linalg.norm(p3)
    
    if p_magnitude == 0:
        return 0.0  # undefined, but set to 0 by convention
    else:
        theta = np.arccos(p3[2] / p_magnitude)
        return theta

def thrust_axis(p3: np.ndarray,
                eps: float = 1e-12,
                include_met: bool = False) -> tuple[np.ndarray, float]:
    """
    Compute the thrust axis  a_T  (unnormalised) and the thrust value  T.

    Parameters
    ----------
    p3 : (n,3) ndarray
        Cartesian 3-momenta of the n visible particles.
    eps : float
        Numerical tolerance used to skip (almost) parallel pairs.
    include_met : bool, default False
        If True, add "-Σ p_i" as an extra  3-vector before maximising.
        This reproduces the Herwig 'doMET' option.

    Returns
    -------
    axis : (3,) ndarray
        a_T  multiplied by the thrust value (so |axis| == T).
    T : float
        Thrust ( 0.5 ≤ T ≤ 1 ).
    """
    # ------------------------------------------------------------------
    # Optionally append the missing-momentum vector
    # ------------------------------------------------------------------
    if include_met:
        met = -p3.sum(axis=0)          # -Σ p
        p3  = np.vstack((p3, met[None, :]))

    n = len(p3)
    if n < 2:
        return np.zeros(3), 0.0

    # ------------------------------------------------------------------
    # Pre-computed global quantities
    # ------------------------------------------------------------------
    p_mag = np.linalg.norm(p3, axis=1)
    p_sum = p_mag.sum()                        # Σ |p|

    # Upper-triangular indices  (i < j)
    i_idx, j_idx = np.triu_indices(n, k=1)

    # Cross-product matrix  cross[i,j] = p_i × p_j
    cross = np.zeros((n, n, 3), dtype=p3.dtype)
    cross[i_idx, j_idx] = np.cross(p3[i_idx], p3[j_idx])
    cross[j_idx, i_idx] = -cross[i_idx, j_idx]            # antisymmetric

    cross_norm = np.linalg.norm(cross, axis=2)            # |p_i×p_j|

    # Dot-product tensor  D[k,i,j] = p_k · (p_i × p_j)
    D = np.einsum('ka,ija->kij', p3, cross)

    # ------------------------------------------------------------------
    # Search for the vector giving the maximal | Σ s_k p_k |
    # ------------------------------------------------------------------
    best_mag2 = 0.0
    best_vec  = np.zeros(3)

    for i, j in combinations(range(n), 2):
        n_ij = cross_norm[i, j]
        if n_ij < eps:                     # p_i  ∥  p_j   → skip pair
            continue

        # Signs for all k w.r.t. the normal to (p_i,p_j)
        sign_all = np.where(D[:, i, j] > 0.0, 1.0, -1.0)

        # Σ s_k p_k  WITHOUT the i,j terms (they will be added explicitly)
        base_sum = (sign_all[:, None] * p3).sum(axis=0) \
                 - sign_all[i] * p3[i] - sign_all[j] * p3[j]

        # Four possible sign combinations for p_i and p_j
        signs = np.array([[-1, -1],
                          [-1,  1],
                          [ 1, -1],
                          [ 1,  1]], dtype=p3.dtype)       # (4,2)

        variants = base_sum + signs @ np.stack((p3[i], p3[j]), axis=0)
        mag2     = np.einsum('ij,ij->i', variants, variants)

        k_max = mag2.argmax()
        if mag2[k_max] > best_mag2:
            best_mag2 = mag2[k_max]
            best_vec  = variants[k_max]

    # ------------------------------------------------------------------
    # Normalise
    # ------------------------------------------------------------------
    best_mag = np.sqrt(best_mag2)
    T = best_mag / p_sum if p_sum > 0 else 0.0
    axis = best_vec / best_mag * T if best_mag > 0 else np.zeros(3)

    return axis, T

# ----------------------------------------------------------------------
# a 4×2 matrix with all ± combinations, created **once** and reused
_SIGN_VARIANTS = np.array([[-1, -1],
                           [-1,  1],
                           [ 1, -1],
                           [ 1,  1]], dtype=np.float32)      # (4,2)
# ----------------------------------------------------------------------

def thrust_belle_std(p3: np.ndarray, eps: float = 1e-12) -> tuple[np.ndarray, float]:
    """
    Based on the C++ ROOT/Belle implementation.
    """
    n = len(p3)
    if n == 0:
        return np.zeros(3), 0.0
    
    # Sum of momentum magnitudes
    p_mags = np.linalg.norm(p3, axis=1)
    sump = p_mags.sum()
    if sump == 0:
        return np.zeros(3), 0.0
    
    best_thrust = 0.0
    best_axis = np.zeros(3)
    
    # Loop over each particle as initial thrust candidate
    for i in range(n):
        rvec = p3[i].copy()
        
        # Ensure z-component is positive (or flip if negative)
        if rvec[2] <= 0.0:
            rvec = -rvec
        
        s = np.linalg.norm(rvec)
        if s == 0.0:
            continue
        rvec = rvec / s  # normalize
        
        # Iterative refinement
        for _ in range(n):  # max n iterations
            rprev = rvec.copy()
            rvec = np.zeros(3)
            
            # Sum momenta aligned with current thrust direction
            for j in range(n):
                qvec = p3[j]
                if np.dot(qvec, rprev) >= 0:
                    rvec += qvec
                else:
                    rvec -= qvec
            
            # Check convergence: if no particle changed sign, break
            converged = True
            for j in range(n):
                qvec = p3[j]
                if np.dot(qvec, rvec) * np.dot(qvec, rprev) < 0:
                    converged = False
                    break
            if converged:
                break
        
        # Calculate thrust value for this candidate
        rvec_mag = np.linalg.norm(rvec)
        if rvec_mag == 0.0:
            continue
            
        ttmp = 0.0
        for j in range(n):
            qvec = p3[j]
            ttmp += abs(np.dot(qvec, rvec))
        
        ttmp /= (sump * rvec_mag)
        rvec = rvec / rvec_mag  # normalize
        
        if ttmp > best_thrust:
            best_thrust = ttmp
            best_axis = rvec.copy()
    
    # Scale axis by thrust value
    best_axis = best_axis * best_thrust
    return best_axis, best_thrust


def thrust_axis_fast(p3: np.ndarray,
                     eps: float = 1e-12,
                     include_met: bool = False
                    ) -> tuple[np.ndarray, float]:
    """
    Same physics as before – just a tighter NumPy implementation.
    """
    if include_met:
        p3 = np.vstack((p3, -p3.sum(axis=0, keepdims=True)))

    n = len(p3)
    if n < 4:
        return thrust_belle_std(p3, eps)

    # ------------------------------------------------------------------  
    # Pre-computation that only needs |p_i| or scalar products
    # ------------------------------------------------------------------  
    p_norm = np.linalg.norm(p3, axis=1)           # |p|
    p_sum  = p_norm.sum()                         # Σ |p|
    #
    # List *all* ordered pairs (i<j) once, then work with
    #  flat arrays instead of doubly-nested loops.
    #
    i_idx, j_idx = np.triu_indices(n, k=1)        # (M,) each
    M = len(i_idx)                                # number of pairs

    # Cross-products for every pair – shape (M,3)
    cross = np.cross(p3[i_idx], p3[j_idx])
    n_ij  = np.linalg.norm(cross, axis=1)         # (M,)

    # Mask out parallel / almost-parallel pairs up-front
    good   = n_ij >= eps
    i_idx, j_idx, cross, n_ij = \
        i_idx[good], j_idx[good], cross[good], n_ij[good]
    if len(i_idx) == 0:
        return np.zeros(3), 0.0
    #
    # ------------------------------------------------------------------  
    #  For all particles k and all surviving pairs ℓ (=0…M′−1)
    #     dots[k,ℓ] =  p_k · (p_i×p_j)[ℓ]
    #  Gives an (n,M′) matrix in one go.
    # ------------------------------------------------------------------  
    dots = p3 @ (cross.T)                          # (n,M′)

    # signs[k,ℓ]  =  ±1  (Eq. in the algorithm)
    signs_all = np.where(dots > 0.0, 1.0, -1.0)    # (n,M′)

    # Σ s_k p_k  for every pair, **without** the i,j terms (3-vector each)
    base_sum = (signs_all[..., None] * p3[:, None, :]).sum(axis=0)  # (M′,3)
    base_sum -= (signs_all[i_idx, np.arange(len(i_idx))][:, None] * p3[i_idx])
    base_sum -= (signs_all[j_idx, np.arange(len(j_idx))][:, None] * p3[j_idx])

    # ------------------------------------------------------------------  
    # Add the four (±p_i ±p_j) possibilities – fully vectorised
    # ------------------------------------------------------------------  
    pij = np.stack((p3[i_idx], p3[j_idx]), axis=1)                # (M′,2,3)
    variants = base_sum[:, None, :] + (_SIGN_VARIANTS @ pij)      # (M′,4,3)
    mag2     = np.einsum('mfc,mfc->mf', variants, variants)       # (M′,4)

    # Pick the best variant for every pair, then the global best
    idx_pair   = np.argmax(mag2, axis=1)                          # (M′,)
    best_pair2 = mag2[np.arange(len(mag2)), idx_pair]             # |Σ s p|²
    k_best     = np.argmax(best_pair2)
    best_vec   = variants[k_best, idx_pair[k_best]]

    # ------------------------------------------------------------------  
    # Normalise
    # ------------------------------------------------------------------  
    best_mag = np.sqrt(best_pair2[k_best])
    T    = best_mag / p_sum if p_sum > 0 else 0.0
    axis = best_vec / best_mag * T if best_mag > 0 else np.zeros(3)

    return axis, T

def thrust_axis_fast_optimized(p3: np.ndarray,
                                    eps: float = 1e-12,
                                    include_met: bool = False
                                   ) -> tuple[np.ndarray, float]:
    """
    Ultra-optimized version with minimal memory allocations and maximum efficiency.
    """
    if include_met:
        p3 = np.vstack((p3, -p3.sum(axis=0, keepdims=True)))
    
    n = len(p3)
    
    if n < 4:
        return thrust_belle_std(p3, eps)
    
    # Convert to float32 for better performance
    p3 = p3.astype(np.float32, copy=False)
    
    # Pre-compute norms
    p_sum = np.linalg.norm(p3, axis=1).sum()
    
    # Generate indices and compute cross products in one shot
    i_idx, j_idx = np.triu_indices(n, k=1)
    cross = np.cross(p3[i_idx], p3[j_idx])
    
    # Filter parallel pairs
    cross_norms = np.linalg.norm(cross, axis=1)
    good = cross_norms >= eps
    if not np.any(good):
        return np.zeros(3, dtype=np.float32), 0.0
    
    # Apply filter once
    i_idx, j_idx, cross = i_idx[good], j_idx[good], cross[good]
    M = len(i_idx)
    
    # Compute signs
    dots = p3 @ cross.T
    signs = np.where(dots > 0.0, 1.0, -1.0)
    
    # Most critical optimization: efficient base_sum computation
    base_sum = signs.T @ p3  # This is the key optimization!
    
    # Subtract i,j contributions
    base_sum -= signs[i_idx, np.arange(M)][:, None] * p3[i_idx]
    base_sum -= signs[j_idx, np.arange(M)][:, None] * p3[j_idx]
    
    # Compute variants - use original efficient method
    pij = np.stack((p3[i_idx], p3[j_idx]), axis=1)
    variants = base_sum[:, None, :] + (_SIGN_VARIANTS @ pij)
    
    # Find best - use original efficient method
    mag2 = np.einsum('mfc,mfc->mf', variants, variants)
    best_variant = np.argmax(mag2, axis=1)
    best_mag2_per_pair = mag2[np.arange(M), best_variant]
    best_pair = np.argmax(best_mag2_per_pair)
    
    # Extract result
    best_vec = variants[best_pair, best_variant[best_pair]]
    best_mag = np.sqrt(best_mag2_per_pair[best_pair])
    
    # Normalize
    T = best_mag / p_sum if p_sum > 0 else 0.0
    axis = best_vec / best_mag * T if best_mag > 0 else np.zeros(3, dtype=np.float32)
    
    return axis, float(T)


def thrust_theta(axis_vec, thrust_val, fold=True):
    """
    Return polar angle of the thrust axis.

    Parameters
    ----------
    axis_vec   : ndarray shape (3,)   (output[0] of thrust_axis_fast)
    thrust_val : float                (output[1] of thrust_axis_fast)
    fold       : bool                 if True, returns 0–π/2 (use |cosθ|)

    Returns
    -------
    theta  : float  (radians)
    """
    if thrust_val <= 0:
        return np.nan
    cos_th = axis_vec[2] / thrust_val          # z-component of unit axis
    if fold:
        cos_th = abs(cos_th)
    return np.arccos(np.clip(cos_th, -1.0, 1.0))

def heavy_jet_mass(p4: np.ndarray, 
                   thrust_axis: np.ndarray,
                   eps: float = 1e-12) -> float:
    """
    Calculate heavy jet mass using thrust axis.
    
    Heavy jet mass is defined as the greater of the invariant masses 
    of two hemispheres separated by the plane normal to the thrust axis.
    
    Args:
        p4: 4-momentum array of shape (n_particles, 4) where columns are [E, px, py, pz]
        thrust_axis: 3D thrust axis vector (should be normalized)
        eps: Small number for numerical stability
        
    Returns:
        heavy_jet_mass: The larger invariant mass of the two hemispheres
    """
    if len(p4) == 0:
        return 0.0
    
    # Extract 3-momentum and energy
    p3 = p4[:, 1:]  # [px, py, pz]
    E = p4[:, 0]    # energy
    
    # Normalize thrust axis just to be safe
    thrust_norm = np.linalg.norm(thrust_axis)
    if thrust_norm < eps:
        return 0.0
    n_thrust = thrust_axis / thrust_norm
    
    # Project momenta onto thrust axis to determine hemispheres
    projections = p3 @ n_thrust  # p⃗ · n̂_thrust
    
    # Separate into two hemispheres
    hemisphere_pos = projections > 0
    hemisphere_neg = projections <= 0
    
    # Handle edge case where all particles are in one hemisphere
    if not np.any(hemisphere_pos) or not np.any(hemisphere_neg):
        # All particles in one hemisphere - return total invariant mass
        total_p4 = p4.sum(axis=0)
        return np.sqrt(max(0, total_p4[0]**2 - np.sum(total_p4[1:]**2)))
    
    # Calculate invariant mass for each hemisphere
    def invariant_mass(mask):
        if not np.any(mask):
            return 0.0
        hemisphere_p4 = p4[mask].sum(axis=0)  # Sum 4-momenta
        E_sum = hemisphere_p4[0]
        p_sum = hemisphere_p4[1:]
        mass_sq = E_sum**2 - np.sum(p_sum**2)
        return np.sqrt(max(0, mass_sq))  # Protect against numerical errors
    
    mass_pos = invariant_mass(hemisphere_pos)
    mass_neg = invariant_mass(hemisphere_neg)
    
    return max(mass_pos, mass_neg)

def sphericity_nonlinear(px: np.ndarray,
                         py: np.ndarray,
                         pz: np.ndarray):
    """

    It follows the same conventions as the C++ version shown in:
    https://github.com/jingyucms/delphi-nanoaod/blob/main/delphi-analysis/include/DataProcessing/include/sphericityTools.h

    • **all** particles are used (no `pwflag` / charge filtering here);
    • the weight factor is |p|^{(r-2)} = |p|^{0} = 1   →   every particle
      contributes with unit weight;
    • the 3×3 momentum tensor is normalised by  Σ|p|^{r} = Σ|p|² ;
    • eigen-values are returned in **ascending** order  
      λ₁ ≤ λ₂ ≤ λ₃ and the corresponding eigen-vectors  
      **v₁ = column 0**, **v₂ = column 1**, **v₃ = column 2**.

    Parameters
    ----------
    px, py, pz : 1-D `array_like`
        Cartesian momentum components of the event.

    Returns
    -------
    dict with keys
        "l1","l2","l3"            : eigen-values λ₁ ≤ λ₂ ≤ λ₃
        "v1","v2","v3"            : eigen-vectors (numpy 3-vectors)
        "matrix"                  : the normalised 3×3 tensor  S_{ij}
    """
    # --- build the tensor --------------------------------------------------
    px, py, pz = map(np.asarray, (px, py, pz))
    p2   = px*px + py*py + pz*pz
    norm = p2.sum()                       #   Σ |p|²      (r = 2)

    # tensor components
    Sxx = (px*px).sum() / norm
    Syy = (py*py).sum() / norm
    Szz = (pz*pz).sum() / norm
    Sxy = (px*py).sum() / norm
    Sxz = (px*pz).sum() / norm
    Syz = (py*pz).sum() / norm

    S = np.array([[Sxx, Sxy, Sxz],
                  [Sxy, Syy, Syz],
                  [Sxz, Syz, Szz]], dtype=float)

    # --- eigen-decomposition ----------------------------------------------
    # numpy.linalg.eigh already gives λ in ascending order
    eigenvalues, eigenvectors = np.linalg.eigh(S)

    l1, l2, l3 = eigenvalues          #  λ₁ ≤ λ₂ ≤ λ₃
    v1, v2, v3 = (eigenvectors[:,0],  #  smallest-eigen-value direction
                  eigenvectors[:,1],
                  eigenvectors[:,2])  #  largest-eigen-value direction

    return {
        "l1": l1, "l2": l2, "l3": l3,
        "v1": v1, "v2": v2, "v3": v3,
        "matrix": S
    }

def calculate_sphericity(px, py, pz):
    """
    Calculate the sphericity matrix and perform eigen-decomposition to obtain
    eigenvalues and eigenvectors for the standard (nonlinear) sphericity (r=2).
    
    In this version, we assume that all particles are used (i.e. no particle flag filtering),
    and the weight factor is unity (since (p2)^(0) = 1).
    
    Additionally, the function returns the cosine of the polar angle (theta) for the
    first eigenvector. Since eigenvectors are normalized, the cosine of theta is simply
    the z-component of the eigenvector.
    
    Parameters:
      px, py, pz : iterables of floats
          The momentum components for each particle.
    
    Returns:
      A dictionary containing:
        "eigenvalues"  : NumPy array of eigenvalues (sorted in ascending order).
        "eigenvectors" : NumPy array of the corresponding eigenvectors as columns.
        "cos_theta_v1" : Cosine of the polar angle for the first eigenvector.
        "type"         : A string, "nonlinear" (since r is assumed to be 2).
    """
    def p2(px_val, py_val, pz_val):
        return px_val * px_val + py_val * py_val + pz_val * pz_val

    n = len(px)
    # Initialize a 3x3 matrix for the sphericity calculation.
    m = np.zeros((3, 3), dtype=float)
    norm = 0.0

    # Loop over each particle.
    for i in range(n):
        p2_val = p2(px[i], py[i], pz[i])
        factor = 1.0  # For r = 2, the exponent is (2-2)/2 = 0, so factor is 1.
        m[0, 0] += px[i] * px[i] * factor
        m[1, 1] += py[i] * py[i] * factor
        m[2, 2] += pz[i] * pz[i] * factor
        m[1, 0] += px[i] * py[i] * factor
        m[2, 0] += px[i] * pz[i] * factor
        m[1, 2] += py[i] * pz[i] * factor
        # Accumulate norm: for r = 2, it is the sum of p2 values.
        norm += p2_val

    # Normalize the matrix.
    if (norm == 0): print(px, py, pz)
    m = m / norm

    # Symmetrize the matrix explicitly.
    m[0, 1] = m[1, 0]
    m[0, 2] = m[2, 0]
    m[2, 1] = m[1, 2]

    # Compute the eigenvalues and eigenvectors.
    # np.linalg.eigh returns eigenvalues in ascending order.
    eigenvalues, eigenvectors = np.linalg.eigh(m)

    # The first eigenvector is taken as the first column in eigenvectors.
    # Since the eigenvectors are normalized, the cosine theta of the first eigenvector
    # is simply its z-component (the third component).
    cos_theta_v1 = eigenvectors[2, 2]

    return {
        "eigenvalues": eigenvalues,
        "eigenvectors": eigenvectors,  # Each column corresponds to an eigenvector.
        "cos_theta_v1": cos_theta_v1,
        "type": "nonlinear"
    }

def calculate_sphericity_robust(px, py, pz, min_particles=2, regularization=1e-12):
    """
    Calculate the sphericity matrix with robust error handling and numerical stability.
    
    Parameters:
    -----------
    px, py, pz : array-like
        The momentum components for each particle
    min_particles : int
        Minimum number of particles required for calculation
    regularization : float
        Small value added to diagonal for numerical stability
    
    Returns:
    --------
    dict with keys:
        "eigenvalues"  : NumPy array of eigenvalues (sorted in ascending order)
        "eigenvectors" : NumPy array of eigenvectors as columns
        "cos_theta_v1" : Cosine of polar angle for the first eigenvector
        "type"         : "nonlinear"
        "status"       : "success", "failed", or "insufficient_particles"
        "error_msg"    : Error message if calculation failed
    """
    
    # Convert to numpy arrays for safety
    px = np.asarray(px, dtype=np.float64)
    py = np.asarray(py, dtype=np.float64) 
    pz = np.asarray(pz, dtype=np.float64)
    
    n = len(px)
    
    # Check minimum particle requirement
    if n < min_particles:
        return {
            "eigenvalues": np.array([np.nan, np.nan, np.nan]),
            "eigenvectors": np.full((3, 3), np.nan),
            "cos_theta_v1": np.nan,
            "type": "nonlinear",
            "status": "insufficient_particles",
            "error_msg": f"Only {n} particles, need at least {min_particles}"
        }
    
    try:
        # Calculate momentum squared for each particle
        p2_values = px*px + py*py + pz*pz
        
        # Check for zero momentum particles
        valid_particles = p2_values > 1e-15
        if not np.any(valid_particles):
            return {
                "eigenvalues": np.array([np.nan, np.nan, np.nan]),
                "eigenvectors": np.full((3, 3), np.nan),
                "cos_theta_v1": np.nan,
                "type": "nonlinear", 
                "status": "failed",
                "error_msg": "All particles have zero momentum"
            }
        
        # Filter out zero momentum particles
        px_valid = px[valid_particles]
        py_valid = py[valid_particles]
        pz_valid = pz[valid_particles]
        p2_valid = p2_values[valid_particles]
        
        # Build sphericity matrix using vectorized operations
        norm = np.sum(p2_valid)
        
        if norm <= 1e-15:
            return {
                "eigenvalues": np.array([np.nan, np.nan, np.nan]),
                "eigenvectors": np.full((3, 3), np.nan),
                "cos_theta_v1": np.nan,
                "type": "nonlinear",
                "status": "failed", 
                "error_msg": "Total momentum squared is zero"
            }
        
        # Construct sphericity matrix efficiently
        m = np.zeros((3, 3), dtype=np.float64)
        
        # Diagonal elements
        m[0, 0] = np.sum(px_valid * px_valid) / norm
        m[1, 1] = np.sum(py_valid * py_valid) / norm
        m[2, 2] = np.sum(pz_valid * pz_valid) / norm
        
        # Off-diagonal elements
        m[0, 1] = m[1, 0] = np.sum(px_valid * py_valid) / norm
        m[0, 2] = m[2, 0] = np.sum(px_valid * pz_valid) / norm
        m[1, 2] = m[2, 1] = np.sum(py_valid * pz_valid) / norm
        
        # Add regularization to diagonal for numerical stability
        np.fill_diagonal(m, m.diagonal() + regularization)
        
        # Verify matrix properties
        if not np.allclose(m, m.T, rtol=1e-10, atol=1e-15):
            print(f"Warning: Matrix not symmetric, max asymmetry: {np.max(np.abs(m - m.T))}")
            # Force symmetry
            m = (m + m.T) / 2
        
        # Check if matrix is positive semi-definite
        try:
            # Try Cholesky decomposition as a quick PSD check
            np.linalg.cholesky(m + np.eye(3) * 1e-10)
        except np.linalg.LinAlgError:
            print("Warning: Matrix may not be positive semi-definite")
        
        # Primary eigenvalue calculation using numpy
        try:
            eigenvalues, eigenvectors = np.linalg.eigh(m)
        except np.linalg.LinAlgError as e:
            # Fallback to scipy with different algorithm
            print(f"NumPy eigh failed: {e}, trying scipy...")
            try:
                eigenvalues, eigenvectors = scipy_eigh(m, driver='evr')
            except Exception as e2:
                # Final fallback: try with condition number improvement
                print(f"SciPy eigh also failed: {e2}, trying with larger regularization...")
                m_reg = m + np.eye(3) * 1e-8
                try:
                    eigenvalues, eigenvectors = np.linalg.eigh(m_reg)
                except Exception as e3:
                    return {
                        "eigenvalues": np.array([np.nan, np.nan, np.nan]),
                        "eigenvectors": np.full((3, 3), np.nan),
                        "cos_theta_v1": np.nan,
                        "type": "nonlinear",
                        "status": "failed",
                        "error_msg": f"All eigenvalue methods failed: {e3}"
                    }
        
        # Validate eigenvalues
        if np.any(np.isnan(eigenvalues)) or np.any(np.isinf(eigenvalues)):
            return {
                "eigenvalues": eigenvalues,
                "eigenvectors": eigenvectors,
                "cos_theta_v1": np.nan,
                "type": "nonlinear",
                "status": "failed",
                "error_msg": "NaN or infinite eigenvalues"
            }
        
        # Eigenvalues should be non-negative and sum to 1 (trace preservation)
        eigenvalues = np.maximum(eigenvalues, 0)  # Ensure non-negative
        
        # Verify trace preservation (should sum to 1 for normalized sphericity matrix)
        trace_sum = np.sum(eigenvalues)
        expected_trace = np.trace(m)
        
        if abs(trace_sum - expected_trace) > 1e-10:
            print(f"Warning: Trace not preserved. Expected: {expected_trace:.10f}, Got: {trace_sum:.10f}")
        
        # Extract cosine of polar angle for the largest eigenvalue eigenvector
        # (NumPy eigh returns eigenvalues in ascending order, so last one is largest)
        cos_theta_v1 = eigenvectors[2, -1]  # z-component of largest eigenvalue eigenvector
        
        return {
            "eigenvalues": eigenvalues,
            "eigenvectors": eigenvectors,
            "cos_theta_v1": cos_theta_v1,
            "type": "nonlinear",
            "status": "success",
            "error_msg": None
        }
        
    except Exception as e:
        return {
            "eigenvalues": np.array([np.nan, np.nan, np.nan]),
            "eigenvectors": np.full((3, 3), np.nan),
            "cos_theta_v1": np.nan,
            "type": "nonlinear", 
            "status": "failed",
            "error_msg": f"Unexpected error: {str(e)}"
        }

### main sphericity calculator
def calculate_sphericity_with_fallback(px, py, pz):
    """
    Wrapper function that tries the robust calculation and falls back to default values.
    """
    result = calculate_sphericity_robust(px, py, pz)
    
    if result["status"] != "success":
        print(f"Sphericity calculation failed: {result['error_msg']}")
        print(f"Particle count: {len(px)}")
        if len(px) > 0:
            p_mag = np.sqrt(np.array(px)**2 + np.array(py)**2 + np.array(pz)**2)
            print(f"Momentum magnitudes: min={np.min(p_mag):.2e}, max={np.max(p_mag):.2e}, mean={np.mean(p_mag):.2e}")
        
        # Return safe default values
        return {
            "eigenvalues": np.array([0.0, 0.0, 1.0]),  # Linear configuration as default
            "eigenvectors": np.eye(3),  # Identity matrix
            "cos_theta_v1": 0.0,  # Default to perpendicular
            "type": "nonlinear"
        }
    
    return result


## I/O functions ##

def OpenFile(file_in,iodir):
    """  file_in -- Input file name
         iodir   -- 'r' readonly  'r+' read+write """
    try:
        ifile=open(file_in, iodir)
    except:
        print("Could not open file: ",file_in)
        sys.exit(1)
    return ifile

def ReadFilesFromList(infile):

    ifile=OpenFile(infile,'r')
    iline=0

    x = ifile.readline()

    filelist=[]
    while x != "":
        iline+=1
        filename=x.rstrip()

        if len(filename)>0 and filename[0] != "#":
            print(filename)
            filelist.append(filename)

        x = ifile.readline()

    return filelist

## Helpers for matching ##

def match_angular(rec: P4Block, gen: P4Block, r):

    # --- distance matrix -------------------------------------------------
    dot      = rec.p @ gen.p.T
    cos      = np.clip(dot / np.outer(rec.norm, gen.norm), -1., 1.)
    rvals    = np.arccos(cos)                  # 0 … π  (always ≥0)

    # --- cost matrix: ΔR, but “inf” for opposite charge ------------------
    cost = rvals.copy()
    #cost[(rec.q[:, None] * gen.q[None, :]) < 0] = 9999   # forbid q·Q ≤ 0

    # Hungarian assignment
    row, col = linear_sum_assignment(cost)

    # keep only pairs with ΔR < match_r  (charge sign already enforced)
    good = rvals[row, col] < r
    ireco, igen = row[good], col[good]

    # ------------- bookkeeping --------------------------------------------
    imiss  = np.setdiff1d(np.arange(gen.pt.size),  igen, assume_unique=True)
    ifake  = np.setdiff1d(np.arange(rec.pt.size),  ireco, assume_unique=True)

    return ireco, igen, imiss, ifake

def match_correspondence_index(rec: P4Block, gen: P4Block, cat: str):
    """
    Match particles using correspondence indices.
    
    Parameters
    ----------
    rec, gen : P4Block
        Reconstructed and generator particle blocks
    cat : str
        Category ('pho' for photon matching with mass check, 'trk', 'ntr')
    
    Returns
    -------
    ireco, igen : arrays
        Indices of matched reconstructed and generator particles
    imiss : array
        Indices of generator particles without matches (missed)
    ifake : array
        Indices of reconstructed particles without matches (fakes)
    """
    if rec.p.shape[0] == 0 or gen.p.shape[0] == 0:
        return np.array([]), np.array([]), np.arange(gen.p.shape[0]), np.arange(rec.p.shape[0])
    
    # Check if correspondence indices are available
    if rec.idx is None or rec.cspidx is None or gen.idx is None or gen.cspidx is None:
        # Fall back to angular matching if correspondence indices not available
        return match_angular(rec, gen, 0.05)  # Default angular match

    # Create index mapping for gen particles
    gen_idx_to_pos = {idx: i for i, idx in enumerate(gen.idx)}
    rec_idx_to_pos = {idx: i for i, idx in enumerate(rec.idx)}
    
    ireco = []
    igen = []
    
    # Find matches from rec -> gen
    for i, rec_csp in enumerate(rec.cspidx):
        if rec_csp != 0 and rec_csp in gen_idx_to_pos:
            # Check if the gen particle also points back to this rec particle
            gen_pos = gen_idx_to_pos[rec_csp]
            if gen.cspidx[gen_pos] == rec.idx[i]:
                # Additional check for photons: generator particle must have mass = 0
                if cat == 'pho':
                    if abs(gen.m[gen_pos]) < 1e-3:  # Generator particle should be massless
                        ireco.append(i)
                        igen.append(gen_pos)
                elif cat == 'ntr':
                    if abs(gen.m[gen_pos]) >= 1e-3:  # Generator particle should be massless
                        ireco.append(i)
                        igen.append(gen_pos)
                else:
                    ireco.append(i)
                    igen.append(gen_pos)
    
    ireco = np.array(ireco)
    igen = np.array(igen)
    
    # Find missed gen particles (cspidx == 0 or no corresponding rec particle)
    matched_gen = set(igen)
    imiss = np.array([i for i in range(gen.p.shape[0]) if i not in matched_gen])
    
    # Find fake rec particles (cspidx == 0 or no corresponding rec particle)
    matched_rec = set(ireco)
    ifake = np.array([i for i in range(rec.p.shape[0]) if i not in matched_rec])
    
    return ireco, igen, imiss, ifake

## Helpers for multi-d EEC calculation##

def get_flat_bin_index(i_2d, j_2d, ny):
    """
    Convert 2D bin indices to flat bin index using same logic as unfold_2dflattern.py
    i_2d, j_2d: bin indices from ROOT histogram (can be 0 for underflow, nx+1 for overflow)
    ny: total number of bins in y direction (including overflow)
    Returns: 0-based flat index for numpy array
    
    ROOT bins: 0=underflow, 1...nx=normal, nx+1=overflow
    We map this directly: flat_bin = i * ny + j (eij changes fastest)
    """
    return i_2d * ny + j_2d

def calculate_event_eec_histogram(pairs_data, template_hist, nx, ny):
    """
    Calculate single-event EEC histogram eec^(k)
    
    pairs_data: list of tuples [(jacobian_val, eij_val, weight), ...]
    template_hist: ROOT 2D histogram to get binning from
    nx, ny: total bins including overflow/underflow
    
    Returns: event histogram vector eec^(k) (flattened), including overflow bins
    """
    total_bins = nx * ny
    
    # Create event histogram vector eec^(k)
    event_eec = np.zeros(total_bins)
    
    # Fill the event histogram - include ALL bins (no overflow check)
    for jacobian_val, eij_val, weight in pairs_data:
        # Find which bins this entry belongs to (including overflow/underflow)
        i_bin = template_hist.GetXaxis().FindBin(jacobian_val)  # Can be 0 (underflow) or nx+1 (overflow)
        j_bin = template_hist.GetYaxis().FindBin(eij_val)       # Can be 0 (underflow) or ny+1 (overflow)
        
        # Convert to flat bin index (eij changes fastest)
        flat_idx = get_flat_bin_index(i_bin, j_bin, ny)
        if flat_idx < total_bins:  # Safety check
            event_eec[flat_idx] += weight
    
    return event_eec


## Helpers for systematic uncertainty evaluation ##

def calc_scale_factor(pt, mode='none', shift=0.0, scale=1.0):
    """
    Calculate pT-dependent scale factor for momentum bias.
    
    Parameters
    ----------
    pt : float or array
        Transverse momentum value(s)
    mode : str, optional
        Type of scaling to apply:
        - 'none': No pT-dependent scaling (returns 1.0 + shift) * scale
        - 'linear': Linear pT-dependent scaling from 1.0 at pT=30 to 1/0.6 at pT=50
    shift : float, optional
        Constant shift to add (default: 0.0)
        Applied as: (base_scale + shift)
    scale : float, optional
        Constant multiplicative scale (default: 1.0)
        Applied as: result * scale
    
    Returns
    -------
    float or array
        Scale factor(s) to apply to momentum components
    
    Examples
    --------
    # No bias at all
    calc_scale_factor(pt, mode='none', shift=0.0, scale=1.0)  # returns 1.0
    
    # Constant 10% increase for all particles
    calc_scale_factor(pt, mode='none', shift=0.1)  # returns 1.1
    
    # Only pT-dependent linear function
    calc_scale_factor(pt, mode='linear', shift=0.0, scale=1.0)
    
    # Both: linear function AND constant shift (not recommended)
    calc_scale_factor(pt, mode='linear', shift=0.05)
    """
    pt = np.asarray(pt)
    
    if mode == 'none':
        # No pT-dependent scaling - just constant shift and scale
        base_scale = np.ones_like(pt)
    elif mode == 'linear':
        # pT-dependent linear scaling
        #base_scale = np.where(
        #    pt < 30, 
        #    1.0,  
        #    np.where(pt <= 50, pt/30., 1./0.6)
        #)
        bias = 0.001+0.00001*pt**2
        base_scale = 1.0 + bias

        #bias_scale = 0.2 * np.tanh(pt / 100)
        #base_scale = 1.0 + bias_scale   

        #pT_threshold=30.0
        #slope=0.008
        #sharpness=5.0

        #ramp = np.maximum(0, pt - pT_threshold) * slope
        #gate = 0.5 * (1 + np.tanh((pt - pT_threshold) / sharpness))
        #bias = ramp * gate
        #base_scale = 1.0 + bias
    else:
        raise ValueError(f"Unknown mode: {mode}. Choose 'none' or 'linear'")
    
    # Apply shift and scale
    final_scale = (base_scale + shift) * scale
    
    return final_scale

def apply_momentum_bias(px, py, pz, m, mode='linear', shift=0.0, scale=1.0):
    """
    Apply pT-dependent momentum bias to particle 4-vectors.
    
    Parameters
    ----------
    px, py, pz : arrays
        Original 3-momentum components
    m : array
        Particle masses
    mode : str, optional
        Scaling mode ('linear', 'constant', 'none')
    shift : float, optional
        Constant shift to add to scale factor
    scale : float, optional
        Constant multiplicative scale
    
    Returns
    -------
    dict with keys: px_new, py_new, pz_new, e_new, pt_new
    """
    # Calculate original pT
    pt_orig = np.sqrt(px**2 + py**2)
    
    # Get scale factor for each particle
    scale_factor = calc_scale_factor(pt_orig, mode=mode, shift=shift, scale=scale)
    
    # Apply scale to 3-momentum
    px_new = px * scale_factor
    py_new = py * scale_factor
    pz_new = pz * scale_factor
    
    # Calculate new energy with consistent mass
    e_new = np.sqrt(px_new**2 + py_new**2 + pz_new**2 + m**2)
    
    # Calculate new pT
    pt_new = np.sqrt(px_new**2 + py_new**2)
    
    return {
        'px': px_new,
        'py': py_new,
        'pz': pz_new,
        'e': e_new,
        'pt': pt_new
    }

def randomly_drop_particles(px, py, pz, m, drop_fraction=0.0, seed=None):
    """
    Randomly drop a fraction of particles from the list.
    
    Parameters
    ----------
    px, py, pz, m : arrays
        Particle momentum components and masses
    drop_fraction : float
        Fraction of particles to drop (0.0 = keep all, 0.1 = drop 10%)
    seed : int, optional
        Random seed for reproducibility. If None, uses truly random seed.
    
    Returns
    -------
    dict with keys: px, py, pz, m, mask
        Filtered arrays and boolean mask of kept particles
    """
    if drop_fraction <= 0.0 or drop_fraction >= 1.0:
        # No dropping or invalid fraction
        return {
            'px': px,
            'py': py,
            'pz': pz,
            'm': m,
            'mask': np.ones(len(px), dtype=bool)
        }
    
    # If no seed provided, numpy will use a random seed automatically
    # If seed is provided, set it for reproducibility
    if seed is not None:
        np.random.seed(seed)
    
    # Create mask: True = keep particle, False = drop particle
    n_particles = len(px)
    keep_mask = np.random.random(n_particles) > drop_fraction
    
    return {
        'px': px[keep_mask],
        'py': py[keep_mask],
        'pz': pz[keep_mask],
        'm': m[keep_mask],
        'mask': keep_mask
    }

def apply_fake_drop_charged(px, py, pz, m, seed=None):
    """
    Apply fake particle dropping with pT-dependent probability for charged particles.
    Simulates tracking inefficiency that increases with pT.
    
    Dropping probability:
    - pT < 30 GeV: 0% (no dropping)
    - 30 ≤ pT ≤ 50 GeV: linearly increases from 0% to 40%
    - pT > 50 GeV: 40% (constant)
    
    Parameters
    ----------
    px, py, pz, m : arrays
        Particle momentum components and masses
    seed : int, optional
        Random seed for reproducibility. If None, uses random seed.
    
    Returns
    -------
    dict with keys: px, py, pz, m, mask
        Filtered arrays and boolean mask of kept particles
    """
    # Set random seed if provided
    if seed is not None:
        np.random.seed(seed)
    
    # Calculate pT for each particle
    pt = np.sqrt(px**2 + py**2)
    
    # Calculate dropping probability for each particle
    drop_prob = np.where(
        pt < 30,
        0.0,  # No dropping below 30 GeV
        np.where(
            pt <= 50,
            0.4 * (pt - 30) / 20.0,  # Linear from 0% to 40% between 30-50 GeV
            0.4  # Constant 40% above 50 GeV
        )
    )
    
    # Generate random numbers and create keep mask
    random_vals = np.random.random(len(px))
    keep_mask = random_vals > drop_prob  # Keep if random > drop_prob
    
    return {
        'px': px[keep_mask],
        'py': py[keep_mask],
        'pz': pz[keep_mask],
        'm': m[keep_mask],
        'mask': keep_mask
    }

def calc_multiplicity_weight_linear(N_total):
    """
    Calculates a linear event weight pivoted around the point of best agreement.
    """
    N_total = np.asarray(N_total)
    
    # Define the pivot point and the new, gentler slope
    pivot_N = 25.0
    pivot_weight = 1.0
    slope = 0.005  # This is your main tuning parameter
    
    # Point-slope form of a line
    weight = slope * (N_total - pivot_N) + pivot_weight
    
    return np.maximum(0, weight)