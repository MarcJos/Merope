//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

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
    void enlarge(double width);
};

}  // namespace sac_de_billes

#include "../Geometry/Cuboid.ixx"


