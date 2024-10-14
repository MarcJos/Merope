//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../MicroToMesh/SphereMesh.hxx"
#include "../MicroToMesh/CylinderMesh.hxx"
#include "../Mesh/GmshWriter.hxx"
#include "../Voronoi/VoroToMeshGraph.hxx"

namespace merope {
namespace mesh {
namespace generator {
namespace  auxi {

template<class INCLUSION>
mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const INCLUSION& inclusion,
    const Point<3>& L, size_t i) {
    mesh::meshStructure::VoroMesh_UnStructureData<3> polyhedronMeshData{};
    if constexpr (is_same_v<smallShape::ConvexPolyhedronInc<3>, INCLUSION>) {
        polyhedronMeshData = microToMesh::voroTranslater::getRawMeshGraph(inclusion.center, inclusion.getInnerInclusions()[i].faces, L);
    } else if constexpr (is_same_v<smallShape::SphereInc<3>, INCLUSION>) {
        polyhedronMeshData = microToMesh::getRawMeshGraph(inclusion.getInnerInclusions()[i]);
    } else if constexpr (is_same_v<smallShape::CylinderInc<3>, INCLUSION>) {
        polyhedronMeshData = microToMesh::getRawMeshGraph(inclusion.getInnerInclusions()[i]);
    } else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("I can't mesh this type of shape");
    }
    return polyhedronMeshData;
}

}  // namespace auxi
}  // namespace generator
}  // namespace mesh
}  // namespace merope
