//! Copyright : see license.txt
//!
//! \brief A simple grid at 1,2 or 3 dimensions
//

#ifndef _GRID_IXX
#define _GRID_IXX 1


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../VTKinout/VTKRead.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

template<typename T>
void Grid::readVTK(VTKRead& vf) {

    // Read the VTK File
    vector < T > tab(ng);
    vf.readMat(tab);

    // Prepare the list
    typename vector<T>::iterator it = tab.begin();
    unsigned short pmin = *it, pmax = *it;
    for (; tab.end() != it; ++it) {
        if (*it < pmin) pmin = *it;
        if (*it > pmax) pmax = *it;
    }

    // Number of phases
    phases.resize(pmax - pmin + 1);

    unsigned short i, j, k;
    unsigned nyz = (unsigned)nz * ny;
    size_t I;
    for (I = 0, k = 0; k < nz; ++k) {
        size_t cj = k;
        for (j = 0; j < ny; ++j, cj += nz) {
            size_t ci = cj;
            for (i = 0; i < nx; ++i, ++I, ci += nyz) {
                phases[tab[I] - pmin].push_back(ci);
            }
        }
    }
}

} // namespace merope

#endif // _GRID_IXX