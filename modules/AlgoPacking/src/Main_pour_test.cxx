//! Copyright : see license.txt
//!
//! \brief
//
#include "StdHeaders.hxx"

#include "AlgoPacking.hxx"
#include "For_testing.hxx"
#include "SphereManipulator.hxx"
#include "AlgoWP.hxx"
#include "AmbiantSpace.hxx"

using  namespace sac_de_billes;
using namespace std;

int main() {
    algoRSA_aux_test::FillMaxRSA<3>();
    algoRSA_aux_test::FillMaxRSA<2>();
    algoRSA_aux_test::FillMaxRSA_2();
    algoRSA_aux_test::PickOnSphere();
    algoRSA_aux_test::Interface();
    algoRSA_aux_test::NamePhases();
    algoRSA_aux_test::nonRegression_Tore();
    algoRSA_aux_test::nonRegression_Shapes();
    algoRSA_aux_test::variousTests();
    algoRSA_aux_test::MultiTestWP();
    algoRSA_aux_test::MultiTestBool();
    algoRSA_aux_test::volInter();
    algoRSA_aux_test::ExclusionDistance();
    algoRSA_aux_test::TestPerf();
    algoRSA_aux_test::GeomTests();
    return 0;
}
