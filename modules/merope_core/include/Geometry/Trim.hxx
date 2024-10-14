//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../../../AlgoPacking/src/Geometry/HalfSpace.hxx"

namespace sac_de_billes {
namespace geomTools {

//! return a solid obtained by pealing a layer from the initial solid
template<unsigned short DIM, class SOLID>
SOLID trim(const SOLID& solid, double layerWidth);

//! return a solid obtained by enlarging by from the initial solid
template<unsigned short DIM, class SOLID>
SOLID enlarge(const SOLID& solid, double layerWidth);

}  // namespace geomTools
}  // namespace sac_de_billes
#include "../Geometry/Trim.ixx"


