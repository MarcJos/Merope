//! Copyright : see license.txt
//!
//! \brief

#include "Geometry/BasicGeometricOperations.hxx"

namespace sac_de_billes {
namespace geomTools {

double determinant(const Point<2>& pt0, const Point<2>& pt1) {
    return pt0[0] * pt1[1] - pt0[1] * pt1[0];
}

}  // namespace  geomTools
}  // namespace  sac_de_billes


