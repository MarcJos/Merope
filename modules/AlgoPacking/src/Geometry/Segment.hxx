//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef ALGOPACKING_SRC_GEOMETRY_SEGMENT_HXX_
#define ALGOPACKING_SRC_GEOMETRY_SEGMENT_HXX_

#include "../Geometry/GeomTypes.hxx"

namespace sac_de_billes {
using namespace std;

//! implements a segment in R^d
template<unsigned short DIM>
struct Segment: public array<Point<DIM>, 2>{
};

} // namespace sac_de_billes

#endif /* ALGOPACKING_SRC_GEOMETRY_SEGMENT_HXX_ */
