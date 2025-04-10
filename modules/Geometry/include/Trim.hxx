//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "AmbiantSpace.hxx"
#include "HalfSpace.hxx"

namespace merope {
namespace geomTools {

//! return a solid obtained by pealing a layer from the initial solid
template<unsigned short DIM, class SOLID>
SOLID trim(const SOLID& solid, double layerWidth);

//! return a solid obtained by enlarging by from the initial solid
template<unsigned short DIM, class SOLID>
SOLID enlarge(const SOLID& solid, double layerWidth);

}  // namespace geomTools
}  // namespace merope
#include "Trim.ixx"


