#ifndef PHDST_UXLINK_HPP
#define PHDST_UXLINK_HPP

namespace phdst
{
    // The Fortran common this mirrors is
    //   COMMON /UXLINK/ LTEMP(2), LRTOP, LSTOP, LTTOP, LITOP,
    //                   LRTEMP, LRWTMP, LRAWUX, LBKTOP, LORTOP, LRTINT, LDTOP
    // i.e. LTEMP is a 2-element array, not a scalar — we were off by one
    // word before, which made LDTOP read back as whatever was in LRTINT.
    extern "C" struct
    {
        int ltemp[2];
        int lrtop;
        int lstop;
        int lttop;
        int litop;
        int lrtemp;
        int lrwtmp;
        int lrawux;
        int lbktop;
        int lortop;
        int lrtint;
        int ldtop;
    } uxlink_;

    inline int &LTEMP(int i) { return uxlink_.ltemp[i - 1]; }  // 1-based
    inline int &LRTOP = uxlink_.lrtop;
    inline int &LSTOP = uxlink_.lstop;
    inline int &LTTOP = uxlink_.lttop;
    inline int &LITOP = uxlink_.litop;
    inline int &LRTEMP = uxlink_.lrtemp;
    inline int &LRWTMP = uxlink_.lrwtmp;
    inline int &LRAWUX = uxlink_.lrawux;
    inline int &LBKTOP = uxlink_.lbktop;
    inline int &LORTOP = uxlink_.lortop;
    inline int &LRTINT = uxlink_.lrtint;
    inline int &LDTOP = uxlink_.ldtop;
}   

#endif // PHDST_UXLINK_HPP