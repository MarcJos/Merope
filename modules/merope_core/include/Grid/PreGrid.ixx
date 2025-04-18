//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericTools/Loops.hxx"

#include "PreGrid.hxx"


namespace merope {

template<unsigned short DIM>
inline array<double, DIM> vox::get_dx_from(array<size_t, DIM> nbNodes,
    array<double, DIM> L) {
    array<double, DIM> dx;
    for (size_t i = 0; i < DIM; i++) {
        dx[i] = L[i] / nbNodes[i];
    }
    return dx;
}

template<unsigned short DIM>
inline array<double, DIM> vox::get_L_from(array<size_t, DIM> nbNodes,
    array<double, DIM> dx) {
    array<double, DIM> L;
    for (size_t i = 0; i < DIM; i++) {
        L[i] = dx[i] * nbNodes[i];
    }
    return L;
}

// With_dx<DIM>

template<unsigned short DIM>
inline void vox::With_dx<DIM>::setDx(const array<double, DIM>& dx_) {
    this->dx = dx_;
    for (size_t i = 0; i < DIM; i++) {
        this->inverse_dx[i] = 1. / dx[i];
    }
    this->halfDiagVoxel = computeLengthHalfDiagVoxel();
}


// GridParameters


template<unsigned short DIM>
inline void vox::GridParameters<DIM>::set_nbNodes_L(
    const array<size_t, DIM>& nbNodes_, const array<double, DIM>& L_) {
    this->setLength(L_);
    this->setNbNodes(nbNodes_);
    this->setDx(vox::get_dx_from<DIM>(nbNodes_, L_));
    this->setSubGridIndices(create_array<DIM, size_t>(0), nbNodes_);
}

template<unsigned short DIM>
inline vox::GridParameters<DIM>::GridParameters(array<size_t, DIM> nbNodes_, array<double, DIM> dx_) :
    InsideTorus<DIM>(get_L_from<DIM>(nbNodes_, dx_)),
    SubArrayDimensions<DIM>(nbNodes_),
    With_dx<DIM>(dx_) {
    this->set_nbNodes_L(nbNodes_, vox::get_L_from<DIM>(nbNodes_, dx_));
}

//
template<unsigned short DIM>
vox::GridParameters<DIM> vox::create_grid_parameters_N_L(array<size_t, DIM> nbNodes, array<double, DIM> L) {
    return vox::create_grid_parameters_N_L<DIM>(nbNodes, create_array<DIM, size_t>(0), nbNodes, L);
};

template<unsigned short DIM>
vox::GridParameters<DIM> vox::create_grid_parameters_N_L(array<size_t, DIM> nbNodes, array<size_t, DIM> nMin, array<size_t, DIM> nMax, array<double, DIM> L) {
    GridParameters<DIM> res(nbNodes, L);
    res.set_nbNodes_L(nbNodes, L);
    res.setSubGridIndices(nMin, nMax);
    return res;
};

// homogenization

template<unsigned short DIM>
inline array<array<long, 2>, DIM> vox::auxi::computeGridLimits(
    const Cuboid<DIM>& cuboid, const auto& grid) {
    array<array<long, 2>, DIM> limits;
    const auto& inverse_dx = grid.getInverseDx();
    for (size_t i = 0; i < DIM; i++) {
        limits[i][0] = static_cast<long>(floor(cuboid.x_min[i] * inverse_dx[i]));
        limits[i][1] = static_cast<long>(ceil(cuboid.x_max[i] * inverse_dx[i]) + 1);
    }
    return limits;
}

template<unsigned short DIM>
inline vector<array<array<long, 2>, DIM>> vox::auxi::intersectGridLimits(const array<array<long, 2>, DIM>& gridLimits,
    const GridParameters<DIM>& preSubGrid) {
    array<vector<array<long, 2>>, DIM> indicesLimits;
    for (size_t i = 0; i < DIM; i++) {
        indicesLimits[i] = vox::auxi::auxi_intersectGridLimits(gridLimits[i], preSubGrid.getNbNodeBigGrid()[i], preSubGrid.getNMin()[i], preSubGrid.getNMax()[i]);
    }
    return cartesianProduct<DIM, array<long, 2>>(indicesLimits);
}

inline vector<array<long, 2>>  vox::auxi::auxi_intersectGridLimits(
    const array<long, 2>& gridLimits, long nbNodes, long nMin,
    long nMax) {
    long ecart = gridLimits[0] - auxi_function::fast_modulo(gridLimits[0], nbNodes);
    nMin += ecart;
    nMax += ecart;
    array<long, 2> currentSegment;
    vector<array<long, 2>> result{};
    while (nMin < gridLimits[1]) {
        currentSegment[0] = max(gridLimits[0], nMin);
        currentSegment[1] = min(gridLimits[1], nMax);
        if (currentSegment[0] < currentSegment[1]) {
            result.push_back(currentSegment);
        }
        //
        nMin += nbNodes;
        nMax += nbNodes;
    }
    return result;
}

}  // namespace merope


