//! Copyright : see license.txt
//!
//! \brief

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Mesh/TestMesh.hxx"
#include "Test/Tests.hxx"
#include "../../AlgoPacking/src/Interface.hxx"
#include "MultiInclusions/LaguerreTess.hxx"
#include "MultiInclusions/MultiInclusions.hxx"
#include "MesoStructure/Structure.hxx"
#include "Voxellation/Voxellizer.hxx"
#include "VTKinout/VTK_adapter.hxx"
#include "Field/CartesianField.hxx"

using namespace merope;
using namespace sac_de_billes;

int main() {
    constexpr unsigned short DIM = 3;
    auto seed = 0;  // Random number generations seed for seperating spheres positionning
    double l3D = 1;  // RVE dimensions
    array<double, DIM> L = { l3D, l3D, l3D };
    vector<array<double, 2>> desiredRPhi = { {0.04,  0.1} };
    auto theSpheres = algoSpheres::throwSpheres<3>(algoSpheres::TypeAlgo::RSA,
        AmbiantSpace::NameShape::Tore, L, seed, desiredRPhi, vector<PhaseType>{1});

    LaguerreTess<DIM> lag(L, theSpheres);
    lag.setAspRatio(array<double, DIM> { 2, 1., 1.});
    MultiInclusions<DIM> mi{};
    mi.setInclusions(lag);
    auto identifiers = mi.getAllIdentifiers();
    mi.changePhase(identifiers, identifiers);

    Structure<DIM> structure(mi);


    // texturing
    // building the map
    mt19937 engine(seed);
    std::map<Identifier, Point<DIM>> identifier_to_orientation{};
    for (size_t id : identifiers) {
        identifier_to_orientation[id] = sac_de_billes::randomShooter::pickOnSphere<DIM>(engine);
    }
    // texturer
    double period = 0.005;
    auto texturer = [period, &identifier_to_orientation](Point<DIM> x, PhaseType phase) {
        return std::sin(geomTools::prodScal<DIM>(x, identifier_to_orientation[phase]) / period);
        };

    //! voxellize basic structure
    size_t n = 128;
    array<size_t, DIM> nbNodes = { n, n, n };

    auto grid_parameters = vox::create_grid_parameters_N_L<DIM>(nbNodes, L);
    auto phaseGrid = vox::voxellizer::transformStructIntoGrid<DIM, vox::VoxelRule::Center>(structure, grid_parameters);
    vtk_adapter::printVTK<unsigned short>(phaseGrid, "Laguerre.vtk");

    //! apply texture
    cerr << phaseGrid[0] << endl;
    vox::CartesianGrid<DIM, double> fieldTextured =
        vox::voxellizer::apply_texture<DIM>(phaseGrid, texturer);
    cerr << fieldTextured[0] << endl;
    vtk_adapter::printVTK<double>(fieldTextured, "textured.vtk");

    return EXIT_SUCCESS;
}
