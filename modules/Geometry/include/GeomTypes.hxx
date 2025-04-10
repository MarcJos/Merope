//! Copyright : see license.txt
//!
//! \briefBasic geometric types
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "GeomConstants.hxx"

namespace merope {

//!
using Identifier = long;

// auxiliary
using PhaseType = long;

// Geometric types
// 1D objects
template<unsigned short DIM>
using DiscPoint = std::array<long, DIM>;

template<unsigned short DIM>
using Point = std::array<double, DIM>;

// checks whether the type T can be understood as a point
template<class T, unsigned short DIM>
constexpr bool is_Point = is_base_of_v<Point<DIM>, T> or is_convertible_v<Point<DIM>, T>;

template<unsigned short DIM>
class RenormPoint;


// 2D objects
template<unsigned short DIM>
class Segment;

template<unsigned short DIM>
class HalfSpace;

// 3D objects
template<unsigned short DIM>
class ConvexPolyhedron;

template<unsigned short DIM>
class SpheroPolyhedron;

template<unsigned short DIM>
class Sphere;

template<unsigned short DIM>
class Ellipse;

template<unsigned short DIM, typename>
class Cylinder;

template<unsigned short DIM>
class Cuboid;

}  // namespace merope


