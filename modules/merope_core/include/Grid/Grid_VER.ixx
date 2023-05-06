//! Copyright : see license.txt
//!
//! \brief 
//!
#ifndef GRID_GRID_VER_IXX_
#define GRID_GRID_VER_IXX_


#include "../../../AlgoPacking/src/Loops.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"
#include "../VTKinout/VTKRead.hxx"
#include "../VTKinout/VTKStream.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {

///Grid_VER
template<unsigned short DIM>
inline void Grid_VER::fromGridPhase(const vox::GridPhase<DIM>& gridPhase, vector<double>& coefficients) {
    this->phases = vector<NodesList>{ };
    this->set_NbPhase(*max_element(gridPhase.begin(), gridPhase.end()) + 1);
    this->set_L<DIM>(vox::get_L_from<DIM>(gridPhase.getNbNodeSubGrid(), gridPhase.getDx()));
    this->set_Nb<DIM>(gridPhase.getNbNodeSubGrid());
    // copies the phase table
    for (size_t i = 0; i < this->getNg(); i++) {
        this->set_Voxel(i, gridPhase[i]);
    }
    this->removeUnusedPhase(coefficients);
}

template<unsigned short DIM>
size_t Grid_VER::get_linear_index(const array<size_t, DIM>& ijk) const {
    static_assert(DIM == 2 or DIM == 3);
    if constexpr (DIM == 3) {
        return vox::auxi::get_linear_index<DIM>(ijk, nbNodes);
    }
    else if constexpr (DIM == 2) {
        return vox::auxi::get_linear_index<DIM>(ijk, array<size_t, DIM>{nbNodes[0], nbNodes[1]});
    }
}

/// homogenization


template<unsigned short DIM>
void Grid_VER::set_L(array<double, DIM> L_) {
    d = DIM;
    lx = L_[0];
    if (DIM == 2 or DIM == 3) {
        ly = L_[1];
    }
    if (DIM == 3) {
        lz = L_[2];
    }
    initCommon();
}

template<unsigned short DIM>
void Grid_VER::set_Nb(array<size_t, DIM> nb) {
    nx = nb[0];
    if (DIM == 2 or DIM == 3) {
        ny = nb[1];
    }
    else {
        ny = 1;
    }
    if (DIM == 3) {
        nz = nb[2];
    }
    else {
        nz = 1;
    }
    initCommon();
}

template<unsigned short DIM>
array<size_t, DIM> Grid_VER::get_coord_index(size_t i) const {
    array<size_t, 3> res;
    res[0] = i / (nbNodes[1] * nbNodes[2]);
    i = i % (nbNodes[1] * nbNodes[2]);
    res[1] = i / (nbNodes[2]);
    i = i % nbNodes[2];
    res[2] = i;
    if constexpr (DIM == 3) {
        return res;
    }
    else if constexpr (DIM == 2) {
        return array<size_t, 2>{res[0], res[1]};
    }
}

template<unsigned short DIM>
inline void Grid_VER::symmetrize(array<size_t, DIM> symDirections) {
    if (DIM != 3) {
        throw runtime_error("Only for DIM=3");
    }
    array<size_t, DIM> nVox, new_nVox, nbCopies;
    array<double, DIM> newL;
    for (size_t i = 0; i < DIM; i++) {
        nVox[i] = this->nbNodes[i];
        nbCopies[i] = auxi_function::puissance<size_t>(2, symDirections[i]);
        new_nVox[i] = nbCopies[i] * this->nbNodes[i];
        newL[i] = nbCopies[i] * this->L[i];
    }
    array<size_t, DIM> NProjec;
    for (auto& list : this->phases) {
        auto originalList = list;
        list = { };
        for (const size_t point : originalList) {
            for (const auto& ijk : getAllIndices_from0<DIM, size_t>(nbCopies)) {
                NProjec = this->get_coord_index<DIM>(point);
                for (size_t i = 0; i < DIM; i++) {
                    NProjec[i] = vox::aux::symmetriseAuxi(NProjec[i], ijk[i], nVox[i]);
                }
                list.push_back(vox::auxi::get_linear_index<DIM>(NProjec, new_nVox));
            }
        }
    }
    this->set_L<DIM>(newL);
    this->set_Nb<DIM>(new_nVox);
}


template<unsigned short DIM>
void vox::symmetrize(Grid_VER& myGrid, array<size_t, DIM> symDirections) {
    assert(myGrid.getDim() == DIM);
    myGrid.symmetrize<DIM>(symDirections);
}

template<unsigned short DIM>
void vox::symmetrize(string inputFileName, string outputFileName, array<size_t, DIM> symDirections) {
    VTKRead inputVTK{ inputFileName.c_str() };
    Grid_VER grid{ inputVTK };
    symmetrize<DIM>(grid, symDirections);
    VTKstream outputVTK{ outputFileName.c_str() };
    grid.toVTKCELL(outputVTK);
}

} // namespace merope

#endif /* GRID_GRID_VER_IXX_ */
