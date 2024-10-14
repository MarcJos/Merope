//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/Geometry/GeomTypes.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
//! todo
//! guarantee output and face_indices have the same order
template<unsigned short DIM>
vector<HalfSpace<DIM>> facesFromVertices(const vector<Point<DIM>>& renormalized_vertices, const vector<vector<long>>& face_indices);

//! todo
template<unsigned short DIM>
vector<Segment<DIM>> edgesFromVertices(const vector<Point<DIM>>& renormalized_vertices, const vector<vector<long>>& face_indices);

}  // namespace merope

#include "AuxiConvexPolyhedron.ixx"


