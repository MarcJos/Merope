//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../Mesh/MeshStructure.hxx"
#include "../Voronoi/VoroInterface.hxx"


namespace merope {
namespace microToMesh {
namespace voroTranslater {

class VoroToMeshGraph : public InsideTorus<3> {
    static constexpr unsigned short DIM = 3;
public:
    //! main constructor
    //! \param L : lengths of the torus
    //! \param centerTesssels : center and weights of the tessels
    VoroToMeshGraph(array<double, DIM> L, const vector<Sphere<DIM>>& centerTessels_, double adimensionnalMergeDistance_ = 1e-5);
    //! return the data related to the voro++
    mesh::meshStructure::VoroMesh_UnStructureData<3> getMeshData() { return rawMeshGraph; }
private:
    //! tolerance for merging points
    double adimensionnalMergeDistance;
    //! internal storage of spheres
    vector<Sphere<DIM>> centerTessels;
    //! internal computations for the gmsh graph
    void auxi_compute();
    //! stores the data for the mesh
    mesh::meshStructure::VoroMesh_UnStructureData<3> rawMeshGraph;
};


namespace auxi {
//! @brief compute the voronoi mesh from analysis of the cells and connectivity of them
//! ensuring each face is correctly shared
//! assumes that the cell identifiers are all different
//! guarantees that solid indices are the same as (cells identifiers + 1)
//! @tparam DIM : dimension = 3
//! @param outputVoroPP : output from Voro++ == tessellation of periodic volume
//! @param adimensionnalMergeDistance : minimal distance to merge 2 points (adimensionnalised vs L)
//! @return : a correct mesh graph
mesh::meshStructure::VoroMesh_UnStructureData<3> computeFromListOfCells(Point<3> L, const vector<merope::voroInterface::SingleCell>& outputVoroPP, double adimensionnalMergeDistance);
}  // namespace  auxi



//! fixme
mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const Point<3>& center, const vector<HalfSpace<3>>& center_to_faces, const Point<3>& L);
//! fixme
mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const voroInterface::SingleCell& singleCell);

}  // namespace voroTranslater
}  // namespace microToMesh
}  // namespace merope


#include "VoroToMeshGraph.ixx"


