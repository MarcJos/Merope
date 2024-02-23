//! Copyright : see license.txt
//!
//! \brief

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Mesh/TestMesh.hxx"
#include "Test/Tests.hxx"

using namespace std;
using namespace sac_de_billes;
using namespace merope;

int main() {
    std::cerr << "---------- Tests::polyCrystalCentroidal(true); -------------" << endl;
    Tests::polyCrystalCentroidal(true);
    std::cerr << "---------- Tests::polyCrystalCentroidal(false); -------------" << endl;
    Tests::polyCrystalCentroidal(false);
    std::cerr << "---------- Tests::polyCrystal_fit_volumes_3D(); -------------" << endl;
    Tests::polyCrystal_fit_volumes_3D("polyCrystal_fit_volumes_3D.vtk");
    std::cerr << "---------- Tests::polyCrystal_fit_volumes_2D(); -------------" << endl;
    Tests::polyCrystal_fit_volumes_2D("polyCrystal_fit_volumes_2D.vtk");
    return EXIT_SUCCESS;
}
