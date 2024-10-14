//! Copyright : see license.txt
//!
//! \brief

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Voronoi/VoroToMeshGraph.hxx"
#include "../../AlgoPacking/src/Geometry/GeomTools_1.hxx"
#include "Mesh/GeoObjects.hxx"
#include "Mesh/MeshStructure.hxx"

#include "container.hh"
#include "pre_container.hh"


#include "MeropeNamespace.hxx"


namespace merope {

namespace microToMesh {
namespace voroTranslater {

using namespace mesh::geoObjects;

mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(array<double, 3> L, const vector<Sphere<3>>& centerTessels) {
    VoroToMeshGraph voroToMesh(L, centerTessels);
    return voroToMesh.getMeshData();
}

mesh::meshStructure::VoroMesh_UnStructureData<3> VoroToMeshGraph::getMeshData() {
    this->auxi_compute();
    return rawMeshGraph;
}

void VoroToMeshGraph::auxi_compute() {
    rawMeshGraph.L = this->getL();
    rawMeshGraph.reset();
    auto listCells = this->getSingleCells();
    size_t nbCells = listCells.size();
    size_t j = 0;
    for (const auto& singleCell : listCells) {
        auto simplifiedCell = mesh::meshStructure::VoroMesh_UnStructureData<DIM>(getRawMeshGraph(singleCell));
        rawMeshGraph.append_with_shift(simplifiedCell);
        cerr << "Built cell" << j << " / " << nbCells << endl;
        j++;
    }
}

mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const voroInterface::SingleCell& singleCell) {
    mesh::meshStructure::VoroMesh_UnStructureData<3> result{ create_array<3>(0.), {}, {}, {}, {}, {}, {}, {} };
    for (size_t i = 0; i < singleCell.vertices.size(); i++) {
        const auto& pt = singleCell.vertices[i];
        result.vecPoint.push_back(mesh::geoObjects::GeoPoint<3>(i, pt));
    }
    long k = 1;
    vector<Identifier> leaves{};
    vector<Identifier> surfaceLoopLeaves{};
    for (size_t index_surf = 0; index_surf < singleCell.faceVertices.size(); index_surf++) {
        const auto& vertices = singleCell.faceVertices[index_surf];
        vector<Identifier> curveLoopLeaves{};
        for (size_t j = 0; j < vertices.size(); j++) {
            size_t j_ = (j < vertices.size() - 1) ? j + 1 : 0;
            //
            leaves = { vertices[j], vertices[j_] };
            result.vecEdge.push_back(mesh::geoObjects::Edge(k, leaves));
            curveLoopLeaves.push_back(k);
            k++;
            //
        }
        result.vecCurveLoop.push_back(mesh::geoObjects::CurveLoop(index_surf + 1, curveLoopLeaves));
        //
        leaves = { static_cast<Identifier>(index_surf + 1) };
        result.vecSurface.push_back(mesh::geoObjects::Surface(index_surf + 1, leaves));
        surfaceLoopLeaves.push_back(index_surf + 1);
        //
    }
    result.vecSurfaceLoop.push_back(mesh::geoObjects::SurfaceLoop(1, surfaceLoopLeaves));
    //
    leaves = { 1 };
    result.vecSolid.push_back(mesh::geoObjects::Solid(1, leaves));
    //
    return result;
}

mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const Point<3>& center, const vector<HalfSpace<3> >& center_to_faces, const Point<3>& L) {
    voro::voronoicell_neighbor voroCell{};
    voroCell.init(-L[0] * 0.5, L[0] * 0.5, -L[1] * 0.5, L[1] * 0.5, -L[2] * 0.5, L[2] * 0.5);
    // build the cell
    for (const auto& face : center_to_faces) {
        voroCell.plane(face.vec()[0], face.vec()[1], face.vec()[2], 2 * face.c());
    }
    //
    auto singleCell = voroInterface::SingleCell(1, center, &voroCell);
    return getRawMeshGraph(singleCell);
}

}  // namespace voroTranslater
}  // namespace microToMesh
}  // namespace merope

