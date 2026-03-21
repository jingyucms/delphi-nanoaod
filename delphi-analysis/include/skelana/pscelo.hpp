#ifndef SKELANA_PSCELO_HPP
#define SKELANA_PSCELO_HPP

#include "skelana/mtrack.hpp"

namespace skelana
{
/* +KEEP,PSCELO. Electron Identification Output information
 *               ( filled by PSFELE from ELEPHANT output )
 *
 * NELOU                - Number of tracks with ELOU info
 * LENELO               - Length of information per electron (34)
 * QELOU(LENELO,MTRACK) - Real array of electron ID output
 * KELOU(LENELO,MTRACK) - Integer array (EQUIVALENCE with QELOU)
 *
 * QELOU( 1,I) - E/p (ECAL energy / track momentum)        EOVERP
 * QELOU( 2,I) - Delta z at shower                          ELDZS
 * QELOU( 3,I) - Delta phi matching (HPC)                   ELDPHI
 * QELOU( 4,I) - Delta phi (FEMC)                           ELDFI
 * QELOU( 5,I) - dE/dx measurement                          ELDEDX
 * QELOU( 6,I) - dE/dx error                                ELEDEDX
 * QELOU( 7,I) - Number of TPC wires for dE/dx              ELNWIR
 * QELOU( 8,I) - Pre-shower E/p                             PREOVP
 * QELOU( 9,I) - Pre-shower shower fraction                 PRSHFI
 * QELOU(10,I) - Pre-shower delta z                         PRDZS
 * QELOU(11,I) - Pre-shower delta phi                       PRDPHI
 * QELOU(12,I) - Pre-shower dE/dx                           PRDEDX
 * QELOU(13,I) - Pre-shower dE/dx expected                  PRDEXP
 * QELOU(14,I) - HPC flag                                   PRHPC
 * QELOU(15,I) - Electron probability                       PROBEL
 * QELOU(16,I) - Number of showers                          PRNSHO
 * QELOU(17,I) - Number of shower/track associations        PRNSOD
 * QELOU(18,I) - RICH gas Cherenkov angle 1                 PRRGA1
 * QELOU(19,I) - RICH gas Cherenkov angle 2                 PRRGA2
 * QELOU(20,I) - Electron flag                              ELEFLG
 * QELOU(21,I) - Electron fit quality                       ELEFIT
 * QELOU(22,I) - Z extrapolation at R=0                     ELZEXR0
 * QELOU(23,I) - Phi extrapolation at R=0                   ELPHIR0
 * QELOU(24,I) - R of conversion vertex                     ELRCON
 * QELOU(25,I) - RICH electron ID                           ELRICH
 * QELOU(26,I) - HPC/FEMC/STIC ECAL flag                   ELHPCF
 * QELOU(27,I) - Pre-shower flag                            ELPREF
 * QELOU(28,I) - Best ECAL energy for electron              EELEC
 * QELOU(29,I) - Inner curvature                            ELCURI
 * QELOU(30,I) - Outer curvature                            ELCURO
 * QELOU(31,I) - Inner chi2                                 ELCH2I
 * QELOU(32,I) - Outer chi2                                 ELCH2O
 * QELOU(33,I) - Rotation angle                             ELROTA
 * QELOU(34,I) - Significance                               ELSIGN
 *
 */

/* +KEEP,PSCELOJJ. Pointers to the El Id output information */
inline const int LENELO  = 34;

inline const int IELOEOP =  1;  // E/p
inline const int IELODZS =  2;  // Delta z at shower
inline const int IELODPH =  3;  // Delta phi (HPC)
inline const int IELODFI =  4;  // Delta phi (FEMC)
inline const int IELODDX =  5;  // dE/dx
inline const int IELOEDD =  6;  // dE/dx error
inline const int IELONWI =  7;  // N wires
inline const int IELOPEP =  8;  // Pre-shower E/p
inline const int IELOPSF =  9;  // Pre-shower shower fraction
inline const int IELOPDZ = 10;  // Pre-shower delta z
inline const int IELOPDP = 11;  // Pre-shower delta phi
inline const int IELOPDD = 12;  // Pre-shower dE/dx
inline const int IELOPID = 13;  // Pre-shower dE/dx expected
inline const int IELOPHP = 14;  // HPC flag
inline const int IELOPEL = 15;  // Electron probability
inline const int IELOPNS = 16;  // N showers
inline const int IELOPNO = 17;  // N shower/track associations
inline const int IELOPG1 = 18;  // RICH gas angle 1
inline const int IELOPG2 = 19;  // RICH gas angle 2
inline const int IELOFLG = 20;  // Electron flag
inline const int IELOEFI = 21;  // Electron fit quality
inline const int IELOZEX = 22;  // Z extrapolation at R=0
inline const int IELOPEX = 23;  // Phi extrapolation at R=0
inline const int IELORCO = 24;  // R of conversion
inline const int IELORIC = 25;  // RICH electron ID
inline const int IELOFES = 26;  // ECAL flag (HPC/FEMC/STIC)
inline const int IELOFPT = 27;  // Pre-shower flag
inline const int IELOBES = 28;  // Best ECAL energy
inline const int IELONIT = 29;  // Inner curvature
inline const int IELONTO = 30;  // Outer curvature
inline const int IELOCIT = 31;  // Inner chi2
inline const int IELOCTO = 32;  // Outer chi2
inline const int IELORAN = 33;  // Rotation angle
inline const int IELOSG2 = 34;  // Significance

extern "C" struct
{
    int nelou;
    int kelou[MTRACK][LENELO];
} pscelo_;

inline int   &NELOU = pscelo_.nelou;
inline int   &KELOU(int i, int j) { return pscelo_.kelou[j - 1][i - 1]; }
inline float &QELOU(int i, int j) { return *reinterpret_cast<float *>(&pscelo_.kelou[j - 1][i - 1]); }

} // namespace skelana

#endif // SKELANA_PSCELO_HPP