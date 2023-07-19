//! Copyright : see license.txt
//!
//! \brief

#ifndef SRC_BASICGEOMETRY_CXX
#define SRC_BASICGEOMETRY_CXX

#include "Geometry/BasicGeometricOperations.hxx"

namespace sac_de_billes {
namespace geomTools {

double determinant(Point<2> pt0, Point<2> pt1) {
    return pt0[0] * pt1[1] - pt0[1] * pt1[0];
}

} // namespace  geomTools
} // namespace  sac_de_billes


#endif // SRC_BASICGEOMETRY_CXX