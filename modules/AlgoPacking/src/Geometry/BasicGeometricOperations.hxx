//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../StdHeaders.hxx"

#include "../Geometry/GeomTypes.hxx"
#include "../Geometry/GeomConstants.hxx"
#include "../AuxiFunctions.hxx"

namespace sac_de_billes {
using namespace std;
namespace geomTools {
//! \briefgeometrical functions

// PROJECTIONS
//! retains only the first coordinates of the oldPoint
template<unsigned short DIM1, unsigned short DIM2>
Point<DIM1> linearProjection(Point<DIM2> oldPoint);
//! projects orthogonally a point x onto a line
template<unsigned short DIM>
Point<DIM> projection(const Segment<DIM>& line, const Point<DIM>& x);

//! does the projection of x inside [0,L)
void projection_periodic_1D(double& x, double L);
void projection_periodic_1D_centered(double& x, const double& L);
//! does the bounced segment [x-r,x+r] inside [0,L).
void bounce_1D(double& x, const double& L, const double& r);

// BASIC LINEAR OPERATIONS
//! \return the scalar product of 2 vectors
template<unsigned short DIM, class C1, class C2>
double prodScal(const C1& v1, const C2& v2);
//! \return the vectorial product of 2 vectors
template<unsigned short DIM, class C1, class C2>
Point<DIM> prodVec(const C1& v1, const C2& v2);
//! \return the Hadamart product
template<unsigned short DIM, class C1, class C2>
Point<DIM> odot(const C1& v1, const C2& v2);
//! @return 2D determinant of 2 points
//! @param pt0, pt1 : ordered 2 2D points
double determinant(const Point<2>& pt0, const Point<2>& pt1);

// NORM AND DISTANCE OPERATIONS
//! \return the squared norm of a vector
template<unsigned short DIM, class T, typename std::enable_if_t<is_Point<T, DIM>, bool> = true >
double normeCarre(const T& v);
//! \return the norm of a vector
template<unsigned short DIM, class T, typename std::enable_if_t<is_Point<T, DIM>, bool> = true >
double norme(const T& v);
//! \return the squared distance between 2 vectors (in the euclidean space)
template<unsigned short DIM, class T, typename std::enable_if_t<is_Point<T, DIM>, bool> = true >
double distanceCarre(const T& v1, const T& v2);
//! renormalizes a vector and returns its previous norm
template<unsigned short DIM>
double renormalize(Point<DIM>& v);
//! \return whether 2 vectors are equal, up to machine error
//! uses relative comparison
template<unsigned short DIM, class T, typename std::enable_if_t<is_Point<T, DIM>, bool> = true >
bool areEqual(const T& v1, const T& v2, double epsilon = geomTools::EPSILON);

}  // namespace geomTools
}  // namespace sac_de_billes

#include "../Geometry/BasicGeometricOperations.ixx"


