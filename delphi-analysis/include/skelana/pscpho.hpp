#ifndef SKELANA_PSCPHO_HPP
#define SKELANA_PSCPHO_HPP

#include "skelana/mtrack.hpp"

namespace skelana
{
    /* +KEEP,PSCPHO.        Photon Identification information in HPC
     *                       ( PA extra-module PHOT (30) )
     *
     *     NPHOT                  - Number of tracks with HPC photons
     *     QPHOT(LENPHO, MTRACK)  - Real    array of PHOT information
     *     KPHOT(LENPHO, MTRACK)  - Integer array of PHOT information
     *
     * Fields (SKELANA v1.05 A.3.14):
     *     QPHOT(1, I) - Energy-weighted shower depth
     *     KPHOT(2, I) - Number of clusters
     *     KPHOT(3, I) - Number of the first layer
     *     KPHOT(4, I) - Number of layers
     *     KPHOT(5, I) - Maximum number of consecutive layers
     *     QPHOT(6, I) - Transverse shower fluctuation measure
     *     QPHOT(7, I) - Longitudinal shower fit (reserved)
     */

    inline const int LENPHO = 7;

    extern "C" struct
    {
        int nphot;
        int kphot[MTRACK][LENPHO];
    } pscpho_;

    inline int   &NPHOT = pscpho_.nphot;
    inline int   &KPHOT(int i, int j) { return pscpho_.kphot[j - 1][i - 1]; }
    inline float &QPHOT(int i, int j) { return *reinterpret_cast<float *>(&pscpho_.kphot[j - 1][i - 1]); }

} // namespace skelana

#endif // SKELANA_PSCPHO_HPP
