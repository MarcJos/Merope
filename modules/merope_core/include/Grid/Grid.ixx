//! Copyright : see license.txt
//!
//! \briefA simple grid at 1,2 or 3 dimensions
//

#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../VTKinout/VTKRead.hxx"

#include "../MeropeNamespace.hxx"
#include "../VTKinout/VTK_adapter.hxx"


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
    auto nyz = nz * ny;
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


template<unsigned short DIM>
void Grid::toVTKCELL_T(VTKstream& fvtk) const {
    vector<unsigned short> ids(ng);
    vtkReorderMaterialIdx(&ids[0]);
    string nameValue = "MaterialId";

    if constexpr (DIM == 2) {
        array<double, DIM> dx = { lx / nx, ly / ny };
        array<size_t, DIM> n = { nx, ny };
        vtk_adapter::printVTK<unsigned short, DIM>(fvtk, n, dx, ids, nameValue);
    } else if constexpr (DIM == 3) {
        array<double, DIM> dx = { lx / nx, ly / ny, lz / nz };
        array<size_t, DIM> n = { nx, ny, nz };
        vtk_adapter::printVTK<unsigned short, DIM>(fvtk, n, dx, ids, nameValue);
    }

}

}  // namespace merope


