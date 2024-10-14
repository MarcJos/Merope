//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Mesh/MeshStructure.hxx"
#include "../Voronoi/VoroInterface.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace microToMesh {
namespace voroTranslater {


class VoroToMeshGraph : public voroInterface::VoroInterface<3> {
    static constexpr unsigned short DIM = 3;
public:
    //! main constructor
    //! \param L : lengths of the torus
    //! \param centerTesssels : center and weights of the tessels
    VoroToMeshGraph(array<double, DIM> L, const vector<Sphere<DIM>>& centerTessels_) : VoroInterface<DIM>(L, centerTessels_), centerTessels(centerTessels_), rawMeshGraph() {}
    //! return the data related to the voro++
    mesh::meshStructure::VoroMesh_UnStructureData<3> getMeshData();
private:
    //! internal storage of spheres
    vector<Sphere<DIM>> centerTessels;
    //! internal computations for the gmsh graph
    void auxi_compute();
    //! stores the data for the mesh
    mesh::meshStructure::VoroMesh_UnStructureData<3> rawMeshGraph;
};

//! fixme
mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(array<double, 3> L, const vector<Sphere<3>>& centerTessels);
//! fixme
mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const Point<3>& center, const vector<HalfSpace<3>>& center_to_faces, const Point<3>& L);
//! fixme
mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const voroInterface::SingleCell& singleCell);

}  // namespace voroTranslater
}  // namespace microToMesh
}  // namespace merope



