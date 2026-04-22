#ifndef SKELANA_PSCVDU_HPP
#define SKELANA_PSCVDU_HPP

#include "skelana/mtrack.hpp"

namespace skelana
{
    /* +KEEP,PSCVDU.                          VD Unassociated hits
     *                                        ( MVDH bank (21) )
     *
     *     NVDHT                - Total number of VD hits
     *     NVDUN                - Number of VD Unassociated hits
     *     NVDUMX               - Max. number of VD  Unass. hits
     *     QVDUN(LENVDU,NVDUMX) - Real    array of VD unass. hits
     *     KVDUN(LENVDU,NVDUMX) - Integer array of VD unass. hits
     *
     *     KVDUN( 1,I) - Module number with the sign of Z
     *     QVDUN( 2,I) - Local X (or Z since 94) coordinate
     *     QVDUN( 3,I) - R coordinate (-R if R-Z is measured)
     *     QVDUN( 4,I) - RPhi (or Z since 94) coordinate
     *     QVDUN( 5,I) - Signal to noise ratio of the hit
     *
     * COMMON /PSCVDU/ NVDHT, NVDUN, KVDUN(LENVDU, NVDUMX)
     * EQUIVALENCE (QVDUN, KVDUN)
     */

    inline const int LENVDU = 5;
    inline const int NVDUMX = 1000;

    extern "C" struct
    {
        int nvdht;
        int nvdun;
        int kvdun[NVDUMX][LENVDU];
    } pscvdu_;

    inline int   &NVDHT = pscvdu_.nvdht;
    inline int   &NVDUN = pscvdu_.nvdun;
    inline int   &KVDUN(int i, int j)   { return pscvdu_.kvdun[j - 1][i - 1]; }
    inline float &QVDUN(int i, int j)   { return *reinterpret_cast<float *>(&pscvdu_.kvdun[j - 1][i - 1]); }

    // Pointer-name constants (from PSCVDUJJ sub-common).
    inline const int IVDUMOD = 1;
    inline const int IVDUXLC = 2;
    inline const int IVDURCO = 3;
    inline const int IVDURPH = 4;
    inline const int IVDUSTN = 5;

} // namespace skelana

#endif // SKELANA_PSCVDU_HPP
