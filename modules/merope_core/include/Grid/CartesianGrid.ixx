//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../../AlgoPacking/include/Voxel.hxx"



namespace merope {
namespace vox {

template<unsigned short DIM, class VOXEL_TYPE>
template<class INT_TYPE>
Point<DIM> CartesianGrid<DIM, VOXEL_TYPE>::getCenterVoxel(const array<INT_TYPE, DIM>& ijk) const {
    Point<DIM> x;
    for (size_t d = 0; d < DIM; d++) {
        x[d] = (ijk[d] + 0.5) * this->getDx()[d];
    }
    return x;
}

}  // namespace vox
}  // namespace merope



