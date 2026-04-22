#ifndef SKELANA_PSCPI0_HPP
#define SKELANA_PSCPI0_HPP

#include "skelana/mtrack.hpp"

namespace skelana
{
    /* +KEEP,PSCPI0.             PI0 Identification information in HPC
     *                            ( PA extra-module PHOT (30) )
     *
     *     NPI0                   - Number of tracks with PI0 info
     *     QPI0(LENPIO, MTRACK)   - Real    array of PI0 information
     *     KPI0(LENPIO, MTRACK)   - Integer array of PI0 information
     *
     * Fields (SKELANA v1.05 A.3.13):
     *     QPI0( 1,I) - Mass of HPCANA tensor fit
     *     QPI0( 2,I) - Rotation angle of the HPCANA tensor fit
     *     QPI0( 3,I) - 1st eigenvalue of the HPCANA tensor fit
     *     QPI0( 4,I) - 2nd eigenvalue of the HPCANA tensor fit
     *     KPI0( 5,I) - Number of connected maxima in the HPCANA
     *     KPI0( 6,I) - Number of expected maxima in the HPCANA
     *     QPI0( 7,I) - Relative energy of the first Gaussian
     *     QPI0( 8,I) - Position of the first Gaussian
     *     QPI0( 9,I) - Width    of the first Gaussian
     *     QPI0(10,I) - Position of the second Gaussian
     *     QPI0(11,I) - Width    of the second Gaussian
     *     QPI0(12,I) - Fitted mass
     *     QPI0(13,I) - 1st eigenvalue
     *     QPI0(14,I) - 2nd eigenvalue
     *     QPI0(15,I) - Rotation angle
     *     QPI0(16,I) - Theta of the shower center
     *     QPI0(17,I) - Phi of the shower center
     *     QPI0(18,I) - Theta of the shower
     *     QPI0(19,I) - Phi of the shower
     *     KPI0(20,I) - Number of OD links
     *     KPI0(21,I) - Number of stray showers
     *     QPI0(22,I) - Chi^2 of the pi0 fit
     *     KPI0(23,I) - Number of cells used for the fit
     *     QPI0(24,I) - Chi^2 of the gamma fit
     *     QPI0(25,I) - Sigma_phi (distribution width)
     *     QPI0(26,I) - Sigma_theta (distribution width)
     */

    inline const int LENPIO = 26;

    extern "C" struct
    {
        int npi0;
        int kpi0[MTRACK][LENPIO];
    } pscpi0_;

    inline int   &NPI0 = pscpi0_.npi0;
    inline int   &KPI0(int i, int j) { return pscpi0_.kpi0[j - 1][i - 1]; }
    inline float &QPI0(int i, int j) { return *reinterpret_cast<float *>(&pscpi0_.kpi0[j - 1][i - 1]); }

} // namespace skelana

#endif // SKELANA_PSCPI0_HPP
