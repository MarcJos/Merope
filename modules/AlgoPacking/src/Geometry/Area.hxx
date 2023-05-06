//! Copyright : see license.txt
//!
//! \brief

#ifndef GEOMETRY_AREA_HXX
#define GEOMETRY_AREA_HXX

#include "../StdHeaders.hxx"

#include "../Geometry/GeomTypes.hxx"
#include "../Geometry/GeomTools_1.hxx"

namespace sac_de_billes {
namespace geomTools {
namespace area {

//! @return the area of a triangle, by the Heron formula
//! @tparam DIM 
//! @param pt0, pt1, pt2  : vertices of the triangle
template<unsigned short DIM>
double triangle(Point<DIM> pt0, Point<DIM> pt1, Point<DIM> pt2);

//! @return : the area of a given polygon
//! @param vertices : the vertices of the polygone, in direct order
//! @warning : not implemented yet
double polygon(vector<Point<2>> vertices);

} // namespace area
} // namespace geomTools
} // namespace  sac_de_billes

#include "Area.ixx"

#endif // GEOMETRY_AREA_HXX
