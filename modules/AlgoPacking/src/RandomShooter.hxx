//! Copyright : see license.txt
//!
//! \brief Contains functions for picking random points
//
#ifndef RANDOMSHOOTER_HXX_
#define RANDOMSHOOTER_HXX_

#include "StdHeaders.hxx"

#include "AmbiantSpace.hxx"

namespace sac_de_billes {
using namespace std;

namespace randomShooter {
//! \return a random point in the cuboid (uniform sampling)
//! \param randGenerator : random generator
//! \param randomReal : uniform distribution on [0,1)
//! \param length : dimensions of the cuboid
template<unsigned short DIM>
Point<DIM> pickInCuboid(mt19937& randGenerator,
        uniform_real_distribution<>& randomReal,
        const Point<DIM>& length);
//! \return a random point on the unit sphere (uniform sampling). Does so by randomly throwing inside the sphere, and then renormalizing.
//! \param randGenerator : random generator
//! \param randomReal : uniform distribution on [0,1)
template<unsigned short DIM>
Point<DIM> pickOnSphere(mt19937& randGenerator,
        uniform_real_distribution<>& randomReal);

} // namespace randomShooter
} // namespace sac_de_billes


#include "RandomShooter.ixx"

#endif /* RANDOMSHOOTER_HXX_ */
