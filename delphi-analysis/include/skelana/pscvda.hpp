#ifndef SKELANA_PSCVDA_HPP
#define SKELANA_PSCVDA_HPP

#include "skelana/mtrack.hpp"

namespace skelana
{
    /* +KEEP,PSCVDA.                          VD Associated hits
     *                                      ( MVDH bank (21) )
     *
     *     NVDAS                     - Number of VD Associated hits
     *     NASHT(MTRACK)             - Number of ass. hits per track
     *     QVDAS(LENVDA,MTRACK,NHIT) - Real    array of VD ass. hits
     *     KVDAS(LENVDA,MTRACK,NHIT) - Integer array of VD ass. hits
     *
     *     KVDAS( 1,I,N) - Module number with the sign of Z
     *     QVDAS( 2,I,N) - Local X (or Z since 94) coordinate
     *     QVDAS( 3,I,N) - R coordinate (-R if R-Z is measured)
     *     QVDAS( 4,I,N) - RPhi (or Z since 94) coordinate
     *     QVDAS( 5,I,N) - Signal to noise ratio of the hit
     *
     * COMMON /PSCVDA/ NVDAS, NASHT(MTRACK), KVDAS(LENVDA, MTRACK, NHIT)
     * EQUIVALENCE (QVDAS, KVDAS)
     */

    inline const int LENVDA = 5;
    inline const int NHIT   = 12;

    extern "C" struct
    {
        int nvdas;
        int nasht[MTRACK];
        int kvdas[NHIT][MTRACK][LENVDA];
    } pscvda_;

    inline int   &NVDAS = pscvda_.nvdas;
    inline int   &NASHT(int i)                 { return pscvda_.nasht[i - 1]; }
    inline int   &KVDAS(int i, int j, int n)   { return pscvda_.kvdas[n - 1][j - 1][i - 1]; }
    inline float &QVDAS(int i, int j, int n)   { return *reinterpret_cast<float *>(&pscvda_.kvdas[n - 1][j - 1][i - 1]); }

    // Pointer-name constants (from PSCVDAJJ sub-common).
    inline const int IVDAMOD = 1;
    inline const int IVDAXLC = 2;
    inline const int IVDARCO = 3;
    inline const int IVDARPH = 4;
    inline const int IVDASTN = 5;

} // namespace skelana

#endif // SKELANA_PSCVDA_HPP
