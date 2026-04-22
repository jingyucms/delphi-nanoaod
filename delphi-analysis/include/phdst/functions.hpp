#ifndef PHDST_FUNCTIONS_HPP
#define PHDST_FUNCTIONS_HPP

#include <string>
#include <cstring>

namespace phdst
{
    extern "C" void phdst_(char *, int *, int *, size_t);
    extern "C" void phset_(char *, int *, size_t);
    extern "C" void phrty_(char *, size_t);
    extern "C" int iphpic_(char *, int *, size_t);
    extern "C" void timed_(float *);
    extern "C" void timex_(float *);

    // LPHPA: Fortran "integer function lphpa(iddp, lin, nump)".
    // Returns the ZEBRA L-address of the sub-bank named iddp under parent bank
    // lin, occurrence number nump (0 = standard for named blocklets). See
    // /cvmfs/delphi.cern.ch/.../src/car/phdstxx.car line 16888 for the full
    // catalogue of accepted names (MAIN / EMNC / HCNC / MUID / ELID / HAID /
    // TRAC / TEOD / TEFA / TEFB / "EMCA.SHOWER" / "EMCA.SHOWER.LAYER", ...).
    extern "C" int lphpa_(const char *iddp, int *lin, int *nump, size_t iddp_len);

    inline void PHDST(const std::string &name, int &&n, int &m)
    {
        char c_name[name.size()];
        size_t len = name.size();
        std::strncpy(c_name, name.c_str(), len);
        phdst_(c_name, &n, &m, len);
    }

    inline void PHSET(const std::string &name, int &&n)
    {
        char c_name[name.size()];
        size_t len = name.size();
        std::strncpy(c_name, name.c_str(), len);
        phset_(c_name, &n, len);
    }

    inline std::string PHRTY()
    {
        char name[] = "1234";
        size_t len = 4;
        phrty_(name, len);

        // strip blanks
        while (len > 0 && name[len - 1] == ' ')
        {
            len--;
        }

        return std::string(name, len);
    }

    inline int IPHPIC(const std::string &name, int &&n)
    {
        char c_name[name.size()];
        size_t len = name.size();
        std::strncpy(c_name, name.c_str(), len);
        return iphpic_(c_name, &n, len);
    }

    inline void TIMED(float &time)
    {
        timed_(&time);
    }

    inline void TIMEX(float &time)
    {
        timex_(&time);
    }

    // C++ wrapper for LPHPA. The Fortran side expects a CHARACTER variable
    // padded to its declared length; the DELPHI convention is 4 characters,
    // but some compound names (e.g. "EMCA.SHOWER") are longer. We pass the
    // raw string and let the caller size it; Fortran LPHPA trims on its end.
    inline int LPHPA(const std::string &name, int lparent, int nump = 0)
    {
        int lin  = lparent;
        int np   = nump;
        return lphpa_(name.c_str(), &lin, &np, name.size());
    }
}

#endif // PHDST_FUNCTIONS_HPP