//! Copyright : see license.txt
//!
//! \brief

#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../MicroToMesh/SphereMesh.hxx"
#include "../MicroToMesh/CylinderMesh.hxx"
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

template<gmsh_writer::MeshMethod meshMethod>
void write_all(std::ostream& f, const mesh::meshStructure::VoroMesh_Periodic_Physical<3>& geoPerStructure) {
    gmsh_writer::auxi::writeDict<meshMethod>("Points", geoPerStructure.dictPoint, f);
    //gmsh_writer::auxi::writeDict("PerPoints", geoPerStructure.dictPerPoint, f);
    gmsh_writer::auxi::writeDict<meshMethod>("Edges", geoPerStructure.dictEdge, f);
    gmsh_writer::auxi::writeDict<meshMethod>("CurveLoops", geoPerStructure.dictCurveLoop, f);
    gmsh_writer::auxi::writeDict<meshMethod>("Surfaces", geoPerStructure.dictSurface, f);
    gmsh_writer::auxi::writeDict<meshMethod>("SurfaceLoop", geoPerStructure.dictSurfaceLoop, f);
    gmsh_writer::auxi::writeDict<meshMethod>("Matrix volume", geoPerStructure.dictSolid, f);
    gmsh_writer::auxi::writeDict<meshMethod>("Periodic surfaces of the enveloppe", geoPerStructure.dictPerSurface, f);
    gmsh_writer::auxi::writeDict<meshMethod>("Physical volumes", geoPerStructure.dictPhysicalVolume, f);
    gmsh_writer::auxi::writeDict<meshMethod>("Physical surfaces", geoPerStructure.dictPhysicalSurface, f);
}


}  // namespace auxi
}  // namespace generator
}  // namespace mesh
}  // namespace merope
