//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace sac_de_billes {
using namespace std;

// define constants
# ifndef M_PI
constexpr double m_PI = 3.1415926535897932384626433832795029L;
#else
constexpr double m_PI = M_PI;
# endif

namespace geomTools {
//! \briefgeometrical functions
//! for comparing and absorbing errors in float operations
constexpr double EPSILON = 1e-6;
constexpr double EPSILON_S = 1e-10;
}  // namespace geomTools
}  // namespace sac_de_billes



