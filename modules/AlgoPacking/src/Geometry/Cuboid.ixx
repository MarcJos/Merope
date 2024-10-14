//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../Geometry/GeomTools_1.hxx"

namespace sac_de_billes {
template<unsigned short DIM>
inline void Cuboid<DIM>::linearTransform(const Point<DIM>& linTransform) {
    linearTransform::point <DIM>(x_max, linTransform);
    linearTransform::point <DIM>(x_min, linTransform);
}

template<unsigned short DIM>
inline void Cuboid<DIM>::enlarge(double width) {
    for (size_t i = 0; i < DIM; i++) {
        x_min[i] -= width;
        x_max[i] += width;
    }
}
}


