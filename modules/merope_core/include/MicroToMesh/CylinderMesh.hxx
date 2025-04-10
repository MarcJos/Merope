//! Copyright : see license.txt
//!
//! \brief

#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/GeomTools_1.hxx"

#include "../Mesh/GeoObjects.hxx"
#include "../Mesh/MeshStructure.hxx"


namespace merope {
namespace microToMesh {

//! fixme
inline mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const Cylinder<3>& cylinder) {
    mesh::meshStructure::VoroMesh_UnStructureData<3> result{ create_array<3>(0.), {}, {}, {}, {}, {}, {}, {}, {} };

    Point<3> basis = cylinder.axis[0];
    // build the local frame
    Point<3> e_3 = RenormPoint<3>(cylinder.axis[1] - cylinder.axis[0]).getPoint();
    Point<3> e_1 = RenormPoint<3>(geomTools::prodVec<3>(Point<3>{ 1., 0, 0 }, e_3)).getPoint();
    if (geomTools::normeCarre<3>(e_1) < 0.1) {
        e_1 = RenormPoint<3>(geomTools::prodVec<3>(Point<3>{ 0, 1, 0 }, e_3)).getPoint();
    }
    Point<3> e_2 = RenormPoint<3>(geomTools::prodVec<3>(e_3, e_1)).getPoint();
    //
    double height = geomTools::norme<3>(cylinder.axis[1] - cylinder.axis[0]);
    double radius = cylinder.radius;

    auto addPoint = [&](Identifier indexLoc, double x, double y, double z) {
        Point<3> point = basis + x * e_1 + y * e_2 + z * e_3;
        result.vecPoint.push_back(mesh::geoObjects::GeoPoint<3>(indexLoc, point));
        };

    // Define the points for the top circle
    addPoint(1, 0, 0, height);
    addPoint(2, radius, 0, height);
    addPoint(3, 0, radius, height);
    addPoint(4, -radius, 0, height);
    addPoint(5, 0, -radius, height);

    // Define the points for the base circle
    addPoint(6, 0, 0, 0);
    addPoint(7, radius, 0, 0);
    addPoint(8, 0, radius, 0);
    addPoint(9, -radius, 0, 0);
    addPoint(10, 0, -radius, 0);

    auto addCircle = [&result](Identifier indexLoc, long i, long j, long k) {
        result.vecEdge.push_back(mesh::geoObjects::Edge(indexLoc, { i, j, k }, mesh::geoObjects::TypeEdge::Circle));
        };
    // Define the lines for the base and top circles
    addCircle(1, 2, 1, 3);
    addCircle(2, 3, 1, 4);
    addCircle(3, 4, 1, 5);
    addCircle(4, 5, 1, 2);

    addCircle(5, 7, 6, 8);
    addCircle(6, 8, 6, 9);
    addCircle(7, 9, 6, 10);
    addCircle(8, 10, 6, 7);


    auto addLine = [&result](Identifier indexLoc, long i, long j) {
        result.vecEdge.push_back(mesh::geoObjects::Edge(indexLoc, { i, j }));
        };
    // Define the lines connecting the base and top circles
    addLine(9, 7, 2);
    addLine(10, 8, 3);
    addLine(11, 9, 4);
    addLine(12, 10, 5);

    auto addLineLoop = [&result](Identifier indexLoc, vector<Identifier> vec) {
        result.vecCurveLoop.push_back(mesh::geoObjects::CurveLoop(indexLoc, vec));
        };
    // Define the surfaces for the cylinder
    addLineLoop(1, { 1, 2, 3, 4 }); // Top circle
    addLineLoop(2, { -8, -7, -6, -5 }); // Top circle
    addLineLoop(3, { 9, -4, -12, 8 }); // Side surface 
    addLineLoop(4, { 10, -1, -9, 5 }); // Side surface 
    addLineLoop(5, { 11, -2, -10, 6 }); // Side surface 
    addLineLoop(6, { 12, -3, -11, 7 }); // Side surface 

    auto addPlaneSurface = [&result](Identifier indexLoc, vector<Identifier> vec) {
        result.vecSurface.push_back(mesh::geoObjects::Surface(indexLoc, vec));
        };

    addPlaneSurface(1, { 1 }); // Top surface
    addPlaneSurface(2, { 2 }); // Base surface

    auto addCurvedSurface = [&result](Identifier indexLoc, vector<Identifier> vec) {
        result.vecSurface.push_back(mesh::geoObjects::Surface(indexLoc, vec, mesh::geoObjects::TypeSurface::Curved));
        };

    addCurvedSurface(3, { 3 }); // Side surface 1
    addCurvedSurface(4, { 4 });
    addCurvedSurface(5, { 5 });
    addCurvedSurface(6, { 6 });


    // Define the volume for the cylinder
    result.vecSurfaceLoop.push_back(mesh::geoObjects::SurfaceLoop(1, { 1, 2, 3, 4, 5, 6 }));
    //
    result.vecSolid.push_back(mesh::geoObjects::Solid(1, { 1 }));
    return result;
}

}  // namespace microToMesh
}  // namespace merope


/* Generalisation of

// Define the radius, height, and mesh size of the cylinder
radius = 1.0;
height = 2.0;
meshSize = 0.1;

// Define the points for the top circle
Point(1) = {0, 0, height, meshSize};
Point(2) = {radius, 0, height, meshSize};
Point(3) = {0, radius, height, meshSize};
Point(4) = {-radius,0, height, meshSize};
Point(5) = {0, -radius, height, meshSize};

// Define the points for the base circle
Point(6) = {0, 0, 0, meshSize};
Point(7) = {radius, 0, 0, meshSize};
Point(8) = {0, radius, 0, meshSize};
Point(9) = {-radius,0, 0, meshSize};
Point(10) = {0, -radius, 0, meshSize};


// Define the lines for the base and top circles
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

Circle(5) = {7, 6, 8};
Circle(6) = {8, 6, 9};
Circle(7) = {9, 6, 10};
Circle(8) = {10, 6, 7};

// Define the lines connecting the base and top circles
Line(9) = {7, 2};
Line(10) = {8, 3};
Line(11) = {9, 4};
Line(12) = {10, 5};


// Define the surfaces for the cylinder
Line Loop(1) = {1, 2, 3, 4}; // Top circle
Line Loop(2) = {-8, -7, -6, -5}; // Top circle
Line Loop(3) = {9, -4, -12, 8}; // Side surface
Line Loop(4) = {10, -1, -9, 5}; // Side surface
Line Loop(5) = {11, -2, -10, 6}; // Side surface
Line Loop(6) = {12, -3, -11, 7}; // Side surface


Plane Surface(1) = {1}; // Top surface
Plane Surface(2) = {2}; // Base surface

Surface(3) = {3}; // Side surface 1
Surface(4) = {4};
Surface(5) = {5};
Surface(6) = {6};


// Define the volume for the cylinder
Surface Loop(1) = {1, 2, 3, 4, 5, 6};
Volume(1) = {1};

// Mesh the volume
Mesh 3;


*/