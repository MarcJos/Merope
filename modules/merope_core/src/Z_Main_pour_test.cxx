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
    Tests::outputLaminate();
    //
    testMesh::test5();
    cerr << "testMesh::test5();" << endl;
    //
    testMesh::test6();
    cerr << "testMesh::test6();" << endl;
    //
    Tests::testFields2();
    cerr << "Tests::testFields2();" << endl;
    testMesh::test3();
    cerr << "testMesh::test3();" << endl;
    testMesh::test4();
    cerr << "testMesh::test4();" << endl;
    testMesh::test2();
    cerr << "testMesh::test2();" << endl;
    testMesh::test1();
    cerr << "testMesh::test1();" << endl;
    ///

    clock_t t_0 = clock();
    cerr << "---4" << endl;
    Tests::extraction();
    cerr << "----5" << endl;

    Tests::polyCrystal11();
    cerr << "-----6" << endl;

    Tests::polyCrystal10();
    cerr << "------7" << endl;

    Tests::polyCrystal0();
    cerr << "--------8" << endl;

    Tests::polyCrystal1();
    cerr << "----------9" << endl;


    Tests::testPerf0();
    cerr << "---------------10" << endl;

    Tests::testFields();
    cerr << "---------------11" << endl;


    cerr << "Time : " << (static_cast<double>(clock() - t_0)) / CLOCKS_PER_SEC << endl;

    return EXIT_SUCCESS;
}
