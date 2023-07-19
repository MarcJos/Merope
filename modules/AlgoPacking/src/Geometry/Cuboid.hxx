//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef ALGOPACKING_SRC_GEOMETRY_CUBOID_HXX_
#define ALGOPACKING_SRC_GEOMETRY_CUBOID_HXX_

#include "../StdHeaders.hxx"

#include "../Geometry/GeomTypes.hxx"

namespace sac_de_billes {
using namespace std;

template<unsigned short DIM>
struct Cuboid {
    //! a cuboid with edges directed by the vectors of the canonical base
    Cuboid(Point<DIM> x_min_, Point<DIM> x_max_) :
        x_min{ x_min_ }, x_max{ x_max_ } {}
    Point<DIM> x_min;
    Point<DIM> x_max;
    //! The space is dilated through (x_1,x_2,x_3) -> (lambda_1 x_1, lambda_2 x_2, lambda_3 x_3)
    void linearTransform(const Point<DIM>& linTransform);
};

} // namespace sac_de_billes

#include "../Geometry/Cuboid.ixx"

#endif /* ALGOPACKING_SRC_GEOMETRY_CUBOID_HXX_ */
