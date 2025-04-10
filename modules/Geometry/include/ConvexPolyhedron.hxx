//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "GeomTypes.hxx"
#include "Cuboid.hxx"

namespace merope {

//! represents a convex polyhedron
template<unsigned short DIM>
class ConvexPolyhedron {
public:
    //! constructor
    ConvexPolyhedron(const Point<DIM>& center__, const vector<HalfSpace<DIM>>& faces__) : center{ center__ }, faces{ faces__ } {}
    //! constructor
    explicit ConvexPolyhedron(const Cuboid<DIM> cuboid);


public:
    //! center of the polyhedron
    //! \warning : origin of the locol coordinates
    Point<DIM> center;
    //! faces of the polyhedron
    //! \warning : in local coordinates (i.e. wrt the center)
    vector<HalfSpace<DIM>> faces;
    //! return whether a point is inside the solid
    bool isInside(const Point<DIM>& point) const;
};

}  // namespace merope

#include "ConvexPolyhedron.ixx"


