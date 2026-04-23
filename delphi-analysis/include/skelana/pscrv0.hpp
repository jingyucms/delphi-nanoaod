#ifndef SKELANA_PSCRV0_HPP
#define SKELANA_PSCRV0_HPP

namespace skelana
{
    /* +KEEP,PSCRV0.                          Reconstructed V0
     *                                        ( V0ID bank (22) )
     *
     *     NRV0                  - Number of reconstructed V0
     *     QRV0(LENRV0, NV0MAX)  - Real    array of V0 information
     *     KRV0(LENRV0, NV0MAX)  - Integer array of V0 information
     *
     * Fields (SKELANA v1.05 A.1.5):
     *     KRV0( 1,I) - Index of the first  particle
     *     KRV0( 2,I) - Index of the second particle
     *     KRV0( 3,I) - Index of the incoming particle (0 if none)
     *     KRV0( 4,I) - Index of the MC incoming particle
     *     KRV0( 5,I) - Tagging flag (bit mask)
     *                    0 - loose, 1 - tight, 2 - K0 bkg, 3 - Lambda bkg,
     *                    22 - K0, 33 - Lambda
     *     QRV0( 6,I) - V0 momentum (GeV/c)
     *     QRV0( 7,I) - Probability of chi2 (PXFVTX fit)
     *     QRV0( 8,I) - X coord. of V0 vertex (cm)
     *     QRV0( 9,I) - Y coord. of V0 vertex (cm)
     *     QRV0(10,I) - Z coord. of V0 vertex (cm)
     *     QRV0(11,I) - xy flight distance normalized to its error
     *     QRV0(12,I) - angle (rad) in xy w.r.t. line to primary vertex
     *     QRV0(13..15,I) - |p+| X/Y/Z after the fit (GeV/c)
     *     QRV0(16..18,I) - |p-| X/Y/Z after the fit (GeV/c)
     *     QRV0(19,I) - suggested V0 mass (negative for antipart.)
     *     QRV0(20,I) - epsilon impact parameter of neutral track
     *     QRV0(21,I) - Z impact parameter of neutral track
     *     QRV0(22,I) - theta of neutral track
     *     QRV0(23,I) - phi of neutral track
     *     QRV0(24..33,I) - weight matrix elements (10 symmetric)
     */

    inline const int NV0MAX = 100;
    inline const int LENRV0 = 33;

    extern "C" struct
    {
        int nrv0;
        int krv0[NV0MAX][LENRV0];
    } pscrv0_;

    inline int   &NRV0 = pscrv0_.nrv0;
    inline int   &KRV0(int i, int j) { return pscrv0_.krv0[j - 1][i - 1]; }
    inline float &QRV0(int i, int j) { return *reinterpret_cast<float *>(&pscrv0_.krv0[j - 1][i - 1]); }

    // Pointer-name constants (from PSCRV0JJ sub-common).
    inline const int IV0ID1 =  1;
    inline const int IV0ID2 =  2;
    inline const int IV0IDI =  3;
    inline const int IV0IDM =  4;
    inline const int IV0FLG =  5;
    inline const int IV0PPP =  6;
    inline const int IV0PRO =  7;
    inline const int IV0XXX =  8;
    inline const int IV0DIS = 11;
    inline const int IV0ANG = 12;
    inline const int IV0PPX = 13;
    inline const int IV0PNX = 16;
    inline const int IV0MAS = 19;

} // namespace skelana

#endif // SKELANA_PSCRV0_HPP
