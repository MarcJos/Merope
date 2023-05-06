//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef ALGOPACKING_SRC_GEOMETRY_CUBOID_IXX_
#define ALGOPACKING_SRC_GEOMETRY_CUBOID_IXX_

#include "../Geometry/GeomTools_1.hxx"

namespace sac_de_billes {
template<unsigned short DIM>
inline void Cuboid<DIM>::linearTransform(const Point<DIM>& linTransform) {
    linearTransform::point <DIM>(x_max, linTransform);
    linearTransform::point <DIM>(x_min, linTransform);
}
}

#endif /* ALGOPACKING_SRC_GEOMETRY_CUBOID_IXX_ */
