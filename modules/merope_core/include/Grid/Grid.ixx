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


template<unsigned short DIM>
void Grid::toVTKCELL_T(VTKstream& fvtk) const {
    // File header
    VTKheaderCELL_T<DIM>(fvtk);
    fvtk << "SCALARS MaterialId unsigned_short" << endl;
    fvtk << "LOOKUP_TABLE default" << endl;
    // Geometry reconstruction
    vector<unsigned short> ids(ng);
    vtkReorderMaterialIdx(&ids[0]);
    // Write phase indices values in file
    fvtk.writeVector(ids, array<size_t, 3>{nx, ny, nz});
}

template<unsigned short DIM>
void Grid::VTKheaderCELL_T(VTKstream& fvtk) const {
    static_assert(DIM == 2 or DIM == 3);
    if constexpr (DIM == 2) {
        const double dx = lx / nx, dy = ly / ny;
        fvtk.STRUCTURED_POINTS(nx, ny, dx, dy);
    } else if constexpr (DIM == 3) {
        const double dx = lx / nx, dy = ly / ny, dz = lz / nz;
        fvtk.STRUCTURED_POINTS(nx, ny, nz, dx, dy, dz);
    }
    // Data type and associated color table
    fvtk.setCELL(ng);
}

} // namespace merope

#endif // _GRID_IXX
