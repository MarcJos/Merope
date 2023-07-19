//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef ALGOPACKING_SRC_GEOMETRY_GEOMCONSTANTS_HXX_
#define ALGOPACKING_SRC_GEOMETRY_GEOMCONSTANTS_HXX_

namespace sac_de_billes {
using namespace std;

// define constants
# ifndef M_PI
constexpr double m_PI = 3.1415926535897932384626433832795029L;
#else
constexpr double m_PI = M_PI;
# endif

namespace geomTools {
//! \brief geometrical functions
//! for comparing and absorbing errors in float operations
constexpr double EPSILON = 1e-6;
constexpr double EPSILON_S = 1e-10;
} // namespace geomTools
} // namespace sac_de_billes


#endif /* ALGOPACKING_SRC_GEOMETRY_GEOMCONSTANTS_HXX_ */
