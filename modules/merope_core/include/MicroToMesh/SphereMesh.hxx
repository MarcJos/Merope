//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "../../AlgoPacking/src/Geometry/GeomTools_1.hxx"
#include "../Mesh/GeoObjects.hxx"
#include "../Mesh/MeshStructure.hxx"
#include "../MeropeNamespace.hxx"


namespace merope {
namespace microToMesh {

//! fixme
inline mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const Sphere<3>& sphere) {
    mesh::meshStructure::VoroMesh_UnStructureData<3> result{ create_array<3>(0.), {}, {}, {}, {}, {}, {}, {} };

    auto easyCoord = [&sphere](const auto& i, const auto& j, const auto& k) {
        return Point<3>({ sphere.center[0] + i * sphere.radius, sphere.center[1] + j * sphere.radius, sphere.center[2] + k * sphere.radius });
        };

    auto writePt = [&result, &easyCoord](Identifier indexLoc, long i, long j, long k) {
        result.vecPoint.push_back(mesh::geoObjects::GeoPoint<3>(indexLoc, easyCoord(i, j, k)));
        };
    //
    writePt(1, 0, 0, 0);
    writePt(2, 1, 0, 0);
    writePt(3, 0, 1, 0);
    writePt(4, 0, 0, 1);
    writePt(5, -1, 0, 0);
    writePt(6, 0, -1, 0);
    writePt(7, 0, 0, -1);
    //

    auto writeCircle = [&result](Identifier indexLoc, long i, long j, long k) {
        result.vecEdge.push_back(mesh::geoObjects::Edge(indexLoc, { i, j, k }, mesh::geoObjects::TypeEdge::Circle));
        };
    //
    writeCircle(1, 2, 1, 3);
    writeCircle(2, 3, 1, 5);
    writeCircle(3, 5, 1, 6);
    writeCircle(4, 6, 1, 2);
    writeCircle(5, 2, 1, 7);
    writeCircle(6, 7, 1, 5);
    writeCircle(7, 5, 1, 4);
    writeCircle(8, 4, 1, 2);
    writeCircle(9, 6, 1, 7);
    writeCircle(10, 7, 1, 3);
    writeCircle(11, 3, 1, 4);
    writeCircle(12, 4, 1, 6);
    //
    auto insertCurveLoop = [&result](Identifier indexLoc, long i, long j, long k) {
        result.vecCurveLoop.push_back(mesh::geoObjects::CurveLoop(indexLoc, { i, j, k }));
        };

    //
    insertCurveLoop(1, 1, 11, 8);
    insertCurveLoop(2, 2, 7, -11);
    insertCurveLoop(3, 3, -12, -7);
    insertCurveLoop(4, 4, -8, 12);
    insertCurveLoop(5, 5, 10, -1);
    insertCurveLoop(6, -2, -10, 6);
    insertCurveLoop(7, -3, -6, -9);
    insertCurveLoop(8, -4, 9, -5);
    //

    vector<long> identifiers{};
    for (long i = 1; i < 9; i++) {
        identifiers.push_back(i);
        result.vecSurface.push_back(mesh::geoObjects::Surface(i, { i }, mesh::geoObjects::TypeSurface::Curved));
    }
    //
    result.vecSurfaceLoop.push_back(mesh::geoObjects::SurfaceLoop(1, identifiers));
    //
    result.vecSolid.push_back(mesh::geoObjects::Solid(1, { 1 }));
    //
    return result;
}

}  // namespace microToMesh
}  // namespace merope