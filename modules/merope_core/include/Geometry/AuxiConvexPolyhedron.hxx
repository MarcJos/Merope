//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MEROPE_CORE_SRC_GEOMETRY_AUXICONVEXPOLYHEDRON_HXX_
#define MEROPE_CORE_SRC_GEOMETRY_AUXICONVEXPOLYHEDRON_HXX_


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

} // namespace merope

#include "AuxiConvexPolyhedron.ixx"

#endif /* MEROPE_CORE_SRC_GEOMETRY_AUXICONVEXPOLYHEDRON_HXX_ */
