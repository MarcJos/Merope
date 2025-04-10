//! Copyright : see license.txt
//!
//! \briefContains functions for picking random points
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "../../Geometry/include/AmbiantSpace.hxx"

namespace sac_de_billes {

namespace randomShooter {
//! \return a random point in the cuboid (uniform sampling)
//! \param randGenerator : random generator
//! \param length : dimensions of the cuboid
template<unsigned short DIM>
Point<DIM> pickInCuboid(mt19937& randGenerator,
        const Point<DIM>& length);
//! \return a random point on the unit sphere (uniform sampling). Does so by randomly throwing inside the sphere, and then renormalizing.
//! \param randGenerator : random generator
template<unsigned short DIM>
Point<DIM> pickOnSphere(mt19937& randGenerator);

}  // namespace randomShooter
}  // namespace sac_de_billes


#include "RandomShooter.ixx"


