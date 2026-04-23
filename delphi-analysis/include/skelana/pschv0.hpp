#ifndef SKELANA_PSCHV0_HPP
#define SKELANA_PSCHV0_HPP

#include "skelana/pscrv0.hpp"  // for NV0MAX

namespace skelana
{
    /* +KEEP,PSCHV0.                    Reconstructed V0 Hypotheses
     *                                        ( V0ID bank (22) )
     *
     *     LENHV0                 - Length of the hypothesis array
     *     QRV0Hn(LENHV0, NV0MAX) - Real    arrays of n-th hypothesis
     *     KRV0Hn(LENHV0, NV0MAX) - Integer arrays of n-th hypothesis
     *
     * n selects the fit hypothesis (SKELANA v1.05 A.1.6):
     *   n=1:  1 VTX fit, no mass constraint
     *   n=2:  2 VTX+DIR fit, no mass constraint
     *   n=3: 11 VTX fit, K0 mass constraint
     *   n=4: 12 VTX+DIR fit, K0 mass constraint
     *   n=5: 21 VTX fit, Lambda mass constraint
     *   n=6: 22 VTX+DIR fit, Lambda mass constraint
     *   n=7: 31 VTX fit, anti-Lambda mass constraint
     *   n=8: 32 VTX+DIR fit, anti-Lambda mass constraint
     *   KRV0Hn(1, I) == -1 if fit failed or low probability
     *
     * Fields (shared across the 8 hypotheses):
     *   KRV0Hn( 1,I) - Kind of fit (see above)
     *   QRV0Hn( 2,I) - Probability of the fit
     *   QRV0Hn( 3,I) - Fitted X coord of the V0 vertex
     *   QRV0Hn( 4,I) - Fitted Y coord of the V0 vertex
     *   QRV0Hn( 5,I) - Fitted Z coord of the V0 vertex
     *   QRV0Hn( 6,I) - Fitted Px of the V0
     *   QRV0Hn( 7,I) - Fitted Py of the V0
     *   QRV0Hn( 8,I) - Fitted Pz of the V0
     *   QRV0Hn( 9,I) - Impact xy
     *   QRV0Hn(10,I) - Error of impact xy
     *   QRV0Hn(11,I) - Impact z
     *   QRV0Hn(12,I) - Error of impact z
     *   QRV0Hn(13,I) - Impact 3D
     *   QRV0Hn(14,I) - Error of impact 3D
     */

    inline const int LENHV0 = 14;

    extern "C" struct
    {
        int krv0h1[NV0MAX][LENHV0];
        int krv0h2[NV0MAX][LENHV0];
        int krv0h3[NV0MAX][LENHV0];
        int krv0h4[NV0MAX][LENHV0];
        int krv0h5[NV0MAX][LENHV0];
        int krv0h6[NV0MAX][LENHV0];
        int krv0h7[NV0MAX][LENHV0];
        int krv0h8[NV0MAX][LENHV0];
    } pschv0_;

    // Hypothesis n is 1..8. Returns a reference into the right slab.
    inline int &KRV0H(int n, int i, int j)
    {
        int (*slabs[8])[LENHV0] = {
            pschv0_.krv0h1, pschv0_.krv0h2, pschv0_.krv0h3, pschv0_.krv0h4,
            pschv0_.krv0h5, pschv0_.krv0h6, pschv0_.krv0h7, pschv0_.krv0h8,
        };
        return slabs[n - 1][j - 1][i - 1];
    }
    inline float &QRV0H(int n, int i, int j)
    {
        return *reinterpret_cast<float *>(&KRV0H(n, i, j));
    }

    // Pointer-name constants (from PSCHV0JJ sub-common).
    inline const int IV0RKI =  1;
    inline const int IV0RPR =  2;
    inline const int IV0RXX =  3;
    inline const int IV0RPX =  6;
    inline const int IV0IXY =  9;
    inline const int IV0DXY = 10;
    inline const int IV0IPZ = 11;
    inline const int IV0DIZ = 12;
    inline const int IV0I3D = 13;
    inline const int IV0D3D = 14;

} // namespace skelana

#endif // SKELANA_PSCHV0_HPP
