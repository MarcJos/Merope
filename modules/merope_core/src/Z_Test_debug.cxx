//! Copyright : see license.txt
//!
//! \brief


#include "Mesh/TestMesh.hxx"
#include "Test/Tests.hxx"
#include "MultiInclusions/SphereInclusions.hxx"
#include "MultiInclusions/MultiInclusions.hxx"
#include "MicroToMesh/MeshGenerator.hxx"
#include "MesoStructure/Structure.hxx"
#include "Voxellation/DynamicVoxellizer.hxx"

using namespace merope;
using namespace sac_de_billes;

int main() {
    Point<3> L = { 1, 1, 1 };
    int nbSpheres = 16;
    double distMin = 0.05;
    int randomSeed = 0;
    auto theSpheres = sac_de_billes::algoSpheres::fillMaxRSA<3>(AmbiantSpace::NameShape::Tore, L, nbSpheres, randomSeed, distMin);
    merope::voroInterface::VoroInterface<3> voro(L, theSpheres, { false, false, false });
    voro.computeSolids();
    return EXIT_SUCCESS;
}