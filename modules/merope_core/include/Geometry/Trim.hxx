//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MEROPE_CORE_SRC_GEOMETRY_TRIM_HXX_
#define MEROPE_CORE_SRC_GEOMETRY_TRIM_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../../../AlgoPacking/src/Geometry/HalfSpace.hxx"

namespace sac_de_billes {
namespace geomTools {

//! return a solid obtained by pealing a layer from the initial solid
template<unsigned short DIM, class SOLID>
SOLID trim(const SOLID& solid, double layerWidth);

} // namespace Corners_of_Cubes
} // namespace sac_de_billes
#include "../Geometry/Trim.ixx"

#endif /* MEROPE_CORE_SRC_GEOMETRY_TRIM_HXX_ */
