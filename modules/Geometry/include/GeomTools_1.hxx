//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "../../GenericTools/AuxiFunctions.hxx"
#include "../../GenericTools/CPP_Functions.hxx"

#include "GeomTypes.hxx"
#include "GeomConstants.hxx"
#include "BasicGeometricOperations.hxx"

namespace merope {

//! \return the number of subcubes when cutting it a Number pieces in any direction
template<unsigned short DIM>
constexpr unsigned short nbSubcubes(unsigned short Number = 2) {
    return auxi_function::puissance<DIM>(Number);
}

namespace geomTools {
//! \briefgeometrical functions

//! compute the volume of a given solid
template<class SOLID>
inline double volume(const SOLID& solid) { return solid.volume(); }
//! \return the volume of many (implicitly non-intersecting) solids
template<class VEC_SOLID>
double volume_all(const VEC_SOLID& s);

//! return whether the point is inside the solid
template<unsigned short DIM, class SOLID>
inline bool isInside(const SOLID& solid, const Point<DIM>& x) { return solid.isInside(x); }
//! \return if inside an intersection of solids
template<unsigned short DIM, class SOLID>
bool isInside_Intersection(const vector<SOLID>& solid, const Point<DIM>& x);
//! \return if inside a union of solids
template<unsigned short DIM, class SOLID>
bool isInside_Union(const vector<SOLID>& solid, const Point<DIM>& x);
//! \return whether the point is on the boundary of a given solid,
//! by testing is point - normal is inside the solid and point + normal is outside it
template<unsigned short DIM, class SOLID>
bool isPointOnBoundary(const SOLID& solid, const Point<DIM>& point, const Point<DIM>& normal, double tolerance) {
    return geomTools::isInside<DIM>(solid, point - tolerance * normal) and not(geomTools::isInside<DIM>(solid, point + tolerance * normal));
}

//! \return signed distance w.r.t to the surface of Solid (<0 = inside, >0 = outside)
//! the point is renormalized wrt the center
template<unsigned short DIM, class SOLID>
double distanceTo(const SOLID& halfSpaces, const Point<DIM>& vector_from_center_to_point);

}  // namespace geomTools

namespace linearTransform {
//! The space is dilated through (x_1,x_2,x_3) -> (lambda_1 x_1, lambda_2 x_2, lambda_3 x_3)
template<unsigned short DIM>
void point(Point<DIM>& vec, const Point<DIM>& linTransform);

template<unsigned short DIM, class C>
double normal(C& normal, const Point<DIM>& linTransform);

template<unsigned short DIM>
double determinant(const Point<DIM>& linTransform);

//! applies the renormalization to Halfplanes
template<unsigned short DIM>
void proceed(HalfSpace<DIM>&, const Point<DIM>& linTransform);
//! applies the renormalization to Spheres
template<unsigned short DIM>
void proceed(Sphere<DIM>&, const Point<DIM>& linTransform);
//! applies the renormalization to Spheres
template<unsigned short DIM>
void proceed(ConvexPolyhedron<DIM>&, const Point<DIM>& linTransform);
//! applies the renormalization to Spheres
template<unsigned short DIM>
void proceed(Ellipse<DIM>&, const Point<DIM>& linTransform);

}  // namespace linearTransform
}  // namespace merope

#include "GeomTools_1.ixx"

