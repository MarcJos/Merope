//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "GeomTools_1.hxx"

namespace merope {

template<unsigned short DIM>
inline merope::ConvexPolyhedron<DIM>::ConvexPolyhedron(
    const Cuboid<DIM> cuboid) : ConvexPolyhedron(0.5 * (cuboid.x_min + cuboid.x_max), {}) {
    Point<DIM> xmin = cuboid.x_min - this->center;
    Point<DIM> xmax = cuboid.x_max - this->center;
    for (size_t i = 0; i < DIM; i++) {
        Point<DIM> outerNormal = create_array<DIM>(0.);
        outerNormal[i] = 1;
        this->faces.push_back(HalfSpace<DIM>(outerNormal, xmax));
        outerNormal[i] = -1;
        this->faces.push_back(HalfSpace<DIM>(outerNormal, xmin));
    }
}

template<unsigned short DIM>
inline bool merope::ConvexPolyhedron<DIM>::isInside(
    const Point<DIM>& point) const {
    return geomTools::isInside_Intersection<DIM>(faces, point - center);
}

}  // namespace merope


