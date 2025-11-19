#ifndef SKELANA_PSCVTX_HPP
#define SKELANA_PSCVTX_HPP

namespace skelana
{
/* +KEEP,PSCVTX.                         Vertex information.
*
*     NVTX                  - Number of reconstr. vertices
*     NVTXMC                - Number of simulated vertices
*     NVTXMX                - Maximum number  of  vertices
*     LENVTX                - Length of vertex information
*     QVTX(LENVTX,2*NVTXMX) - Real    array  of VTX info
*     KVTX(LENVTX,2*NVTXMX) - Integer array  of VTX info
*     LVTX       (2*NVTXMX) - LPV/LSP links for vertices
*
*     KVTX( 1,I) - Index of the first outgoing particle
*     KVTX( 2,I) - Index of the incomming particle
*     KVTX( 3,I) - Nb of outgoing particles (multiplicity)
*     KVTX( 4,I) - Nb of degree of freedom of the vertex fit
*     KVTX( 5,I) - Mass code of the origin particle
*     QVTX( 6,I) - X
*     QVTX( 7,I) - Y    coordinates of the vertex
*     QVTX( 8,I) - Z
*     QVTX( 9,I) - Chi2 of the vertex fit
*     QVTX(10,I) - XX
*     QVTX(11,I) - XY
*     QVTX(12,I) - YY   Error matrix
*     QVTX(13,I) - XZ
*     QVTX(14,I) - YZ
*     QVTX(15,I) - ZZ
*     KVTX(16,I) - Error flag
*     KVTX(17,I) - Vertex status bits :
*                  bit 1 set on if dummy vertex
*                  bit 2 set on if secondary vertex
*                  bit 3 set on if secondary hadronic vertex
*                  bit 4 set on if vertex with simulation data
*
*/

    inline const int LENVTX = 17;
    inline const int NVTXMX = 150;
    extern "C" struct
    {
        int nvtx;
        int nvtxmc;
        int lvtx[2*NVTXMX];
        int kvtx[LENVTX * 2 * NVTXMX];
    } pscvtx_;

    inline int &NVTX = pscvtx_.nvtx;
    inline int &NVTXMC = pscvtx_.nvtxmc;
    inline int& LVTX(int j) {return pscvtx_.lvtx[j-1];}
    inline int& KVTX(int i, int j) {return pscvtx_.kvtx[(i-1) + LENVTX * (j-1)]; }
    inline float& QVTX(int i, int j) {return *reinterpret_cast<float*>(&pscvtx_.kvtx[(i - 1) + LENVTX * (j - 1)]);}
} // namespace skelana

#endif // SKELANA_PSCVTX_HPP
