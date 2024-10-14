//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../Geometry/GeomTypes.hxx"

namespace sac_de_billes {
using namespace std;

//! implements a segment in R^d
template<unsigned short DIM>
struct Segment : public array<Point<DIM>, 2>{
    Point<DIM> middle() const { return 0.5 * ((*this)[0] + (*this)[1]); }
};

}  // namespace sac_de_billes


