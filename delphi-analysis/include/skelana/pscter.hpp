#ifndef SKELANA_PSCTER_HPP
#define SKELANA_PSCTER_HPP

namespace skelana
{
    /* +KEEP,PSCTER.             Unassociated TER's of OD/FCA/FCB
     *                                   ( UTER bank (23) )
     *
     * Track elements (TE ≈ partial-track segments) from three forward
     * sub-detectors that the standard reconstruction did NOT attach to any
     * charged-track object. For particle-flow-style analyses these are the
     * closest thing SKELANA exposes to "hit-level" information (per-cell
     * calorimeter cells are not surfaced by SKELANA — they live in raw SDST
     * banks and would need a separate PHDST-level reader).
     *
     * Same 8-slot layout for OD, FCA, FCB (SKELANA v1.05 A.1.7):
     *     K..( 1,I) - Data descriptor + meas. code  <TER(04)>
     *     Q..( 2,I) - Coordinate 1                  <TER(10)>
     *     Q..( 3,I) - Coordinate 2                  <TER(11)>
     *     Q..( 4,I) - Coordinate 3                  <TER(12)>
     *     Q..( 5,I) - Theta                         <TER(13)>
     *     Q..( 6,I) - Phi                           <TER(14)>
     *     Q..( 7,I) - 1/P or 1/Pt at reference point<TER(15)>
     *     Q..( 8,I) - Length of the track element
     */

    inline const int NUTMAX = 100;
    inline const int LENUTE =   8;

    extern "C" struct
    {
        int nteod;
        int kteod [NUTMAX][LENUTE];
        int ntefca;
        int ktefca[NUTMAX][LENUTE];
        int ntefcb;
        int ktefcb[NUTMAX][LENUTE];
    } pscter_;

    inline int &NTEOD  = pscter_.nteod;
    inline int &NTEFCA = pscter_.ntefca;
    inline int &NTEFCB = pscter_.ntefcb;

    inline int   &KTEOD (int i, int j) { return pscter_.kteod [j - 1][i - 1]; }
    inline float &QTEOD (int i, int j) { return *reinterpret_cast<float *>(&pscter_.kteod [j - 1][i - 1]); }
    inline int   &KTEFCA(int i, int j) { return pscter_.ktefca[j - 1][i - 1]; }
    inline float &QTEFCA(int i, int j) { return *reinterpret_cast<float *>(&pscter_.ktefca[j - 1][i - 1]); }
    inline int   &KTEFCB(int i, int j) { return pscter_.ktefcb[j - 1][i - 1]; }
    inline float &QTEFCB(int i, int j) { return *reinterpret_cast<float *>(&pscter_.ktefcb[j - 1][i - 1]); }

} // namespace skelana

#endif // SKELANA_PSCTER_HPP
