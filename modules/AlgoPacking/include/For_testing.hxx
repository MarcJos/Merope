//! Copyright : see license.txt
//!
//! \briefFor time testing.
//! Comparison with VER in ver/BacASable/Evalue_Temps_VER_TMFFT
#pragma once


#include "../../GenericMerope/StdHeaders.hxx"

#include "../../Geometry/include/AmbiantSpace.hxx"

#include "AlgoPacking.hxx"
#include "AlgoRSA.hxx"
#include "AlgoWP.hxx"
#include "AlgoBool.hxx"
#include "Interface.hxx"
#include "RandomShooter.hxx"
#include "SphereManipulator.hxx"

namespace sac_de_billes {

namespace algoRSA_aux_test {

// for comparison with VER
void inline timeTesting() {
    vector<double> Diameters({ 40., 14., 2. });
    vector<double> Phi2D({ 0.4, 0.0628319, 0.00125664 });
    vector<array<double, 2>> desiredRPhi{ };
    for (size_t i = 0; i < 3; i++) {
        desiredRPhi.push_back(
            array<double, 2> { 0.5 * Diameters[i], Phi2D[i] });
    }
    double lx_0 = 100;
    double ly_0 = 100;

    vector<clock_t> time{ };
    clock_t t_0 = clock();

    for (size_t i = 0; i < 6; i++) {
        double lx = lx_0 * pow(2, i);
        double ly = ly_0 * pow(2, i);
        cout << "lx = " << lx << " ; ly = " << ly << endl;
        AmbiantSpace::Tore<2> T(Point<2> { lx, ly });
        AlgoRSA2D algo(T.L, desiredRPhi, 0, 0, 1, "Tore");

        time.push_back(clock() - t_0);

        cout << "algo           " << "Elapsed CPU time: " << " | nb particles"
            << endl << " 2D, 1st algo :"
            << (static_cast<double>(time[i])) / CLOCKS_PER_SEC << "s |"
            << algo.getSpheres().size() << endl;
        cout << "-------------------" << endl;
    }
}

// for comparison with VER
void inline timeTesting2(string nameFile, double volumeFraction) {
    vector<array<double, 2>> desiredRPhi({
            array<double, 2> { 1, volumeFraction } });
    double lx_0 = 4.;
    double ly_0 = 4.;
    double lz_0 = 4.;

    vector<clock_t> time{ };

    ofstream file(nameFile);

    for (size_t i = 0; i < 6; i++) {
        clock_t t_0 = clock();
        double lx = lx_0 * pow(2, i);
        double ly = ly_0 * pow(2, i);
        double lz = lz_0 * pow(2, i);
        Point<3> L({ lx, ly, lz });
        double exclusionDistance = 0.025;

        AlgoRSA3D algo(L, desiredRPhi, exclusionDistance, 0, 1);

        time.push_back(clock() - t_0);

        cout << lx << " , " << algo.getSpheres().size() << " , "
            << (static_cast<double>(time[i])) / CLOCKS_PER_SEC << endl;

        file << algo.getSpheres().size() << " , "
            << (static_cast<double>(time[i])) / CLOCKS_PER_SEC << endl;
        cout << "-------------------" << endl;
    }
}

inline void nonRegression_Tore() {
    /////////////////////////////////////////////////////////////////////////////////////////////////
    cerr << "----------------------------------------" << endl;
    cerr << "BEGINNING OF NON-REGRESSION TESTS (Tore)" << endl;
    Point<2> TORE_LENGTH_2D({ 64., 64. });
    Point<3> TORE_LENGTH_3D({ 4., 4., 4. });
    double EXCLUSION_DISTANCE = 0;
    double SEED = 2;
    // RESULT
    // 2D, 1st algo :0.911s |26652
    // 2D, 2nd algo :0.418s |26652
    // 3D, 1st algo : 7s | 1770
    // 3D, 2nd algo : 1.9s | 1770
    vector<array<double, 2>> DESIRED_RPHI = { array<double, 2> { 1., 0.25 },
            array<double, 2> { 0.5, 0.1 }, array<double, 2> { 0.25, 0.05 },
            array<double, 2> { 0.125, 1 } };
    AlgoRSA2D* A_2D_1 = new AlgoRSA2D(TORE_LENGTH_2D, DESIRED_RPHI,
        EXCLUSION_DISTANCE, SEED, 1);
    Merope_assert(A_2D_1->getSpheres().size() == 26676,
        "Unexpected number of spheres!");
    /////////////////////////////////////////////////////////////////////////////////////////////////
    delete A_2D_1;
    AlgoRSA2D* A_2D_2 = new AlgoRSA2D(TORE_LENGTH_2D, DESIRED_RPHI,
        EXCLUSION_DISTANCE, SEED, 2);
    Merope_assert(A_2D_2->getSpheres().size() == 26676,
        "Unexpected number of spheres!");
    delete A_2D_2;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    AlgoRSA3D* A_3D_1 = new AlgoRSA3D(TORE_LENGTH_3D, DESIRED_RPHI,
        EXCLUSION_DISTANCE, SEED, 1);
    Merope_assert(A_3D_1->getSpheres().size() == 1760,
        "Unexpected number of spheres!");
    delete A_3D_1;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    AlgoRSA3D* A_3D_2 = new AlgoRSA3D(TORE_LENGTH_3D, DESIRED_RPHI,
        EXCLUSION_DISTANCE, SEED, 2);
    Merope_assert(A_3D_2->getSpheres().size() == 1760,
        "Unexpected number of spheres!");
    delete A_3D_2;
    /////////////////////////////////////////////////////////////////////////////////////////////////
    cerr << "END OF NON-REGRESSION TESTS: SUCCESS!" << endl;
    cerr << "----------------------------------------" << endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////
}

inline void nonRegression_Shapes() {
    /////////////////////////////////////////////////////////////////////////////////////////////////
    cerr << "----------------------------------------" << endl;
    cerr << "BEGINNING OF NON-REGRESSION TESTS (various shapes)" << endl;
    Point<2> ToreLength_2D({ 64., 64. });
    Point<3> ToreLength_3D({ 64., 64., 64. });
    vector<array<double, 2>> desiredRPhi = { array<double, 2> { 1., 1 } };
    unsigned seed = 2;
    double exclusionDistance = 0.0;
    // RESULT
// 	algo           Elapsed CPU time:  | nb particles
//	 2D, 1st algo :0.007s |546
//	 2D, 2nd algo :0.005s | 701
//	 3D, 1st algo : 1.56s | 18249
//	 3D, 2nd algo : 1.478s | 12158
    AlgoRSA2D algo2D_11(ToreLength_2D, desiredRPhi, exclusionDistance, seed, 1,
        "Sphere");
    Merope_assert(algo2D_11.getSpheres().size() == 546,
        "Unexpected number of spheres!");
    AlgoRSA2D algo2D_12(ToreLength_2D, desiredRPhi, exclusionDistance, seed, 2,
        "Sphere");
    Merope_assert(algo2D_12.getSpheres().size() == 546,
        "Unexpected number of spheres!");
    AlgoRSA2D algo2D_21(ToreLength_2D, desiredRPhi, exclusionDistance, seed, 1,
        "Cube");
    Merope_assert(algo2D_21.getSpheres().size() == 701,
        "Unexpected number of spheres!");
    AlgoRSA2D algo2D_22(ToreLength_2D, desiredRPhi, exclusionDistance, seed, 2,
        "Cube");
    Merope_assert(algo2D_22.getSpheres().size() == 701,
        "Unexpected number of spheres!");
    AlgoRSA3D algo3D_11(ToreLength_3D, desiredRPhi, exclusionDistance, seed, 1,
        "Cylinder");
    Merope_assert(algo3D_11.getSpheres().size() == 18286,
        "Unexpected number of spheres!");
    AlgoRSA3D algo3D_12(ToreLength_3D, desiredRPhi, exclusionDistance, seed, 2,
        "Cylinder");
    Merope_assert(algo3D_12.getSpheres().size() == 18286,
        "Unexpected number of spheres!");
    AlgoRSA3D algo3D_21(ToreLength_3D, desiredRPhi, exclusionDistance, seed, 2,
        "Sphere");
    Merope_assert(algo3D_21.getSpheres().size() == 12200,
        "Unexpected number of spheres!");
    AlgoRSA3D algo3D_22(ToreLength_3D, desiredRPhi, exclusionDistance, seed, 2,
        "Sphere");
    Merope_assert(algo3D_22.getSpheres().size() == 12200,
        "Unexpected number of spheres!");
    /////////////////////////////////////////////////////////////////////////////////////////////////
    cerr << "END OF NON-REGRESSION TESTS: SUCCESS!" << endl;
    cerr << "----------------------------------------" << endl;
    /////////////////////////////////////////////////////////////////////////////////////////////////
}

inline void variousTests() {
    cerr << "---------------" << endl;
    cerr << "BEGINNING OF VARIOUS TESTS for RSA" << endl;

    // DATA
    Point<2> ToreLength_2D({ 5., 5. });
    Point<3> ToreLength_3D({ 10., 10., 10. });
    vector<array<double, 2>> desiredRPhi = { array<double, 2> { 1., 0.2 },
            array<double, 2> { 0.5, 0.1 }, array<double, 2> { 0.2, 1 } };
    unsigned seed = 2;
    double exclusionDistance = 0.0;

    clock_t t;
    t = clock();

    vector<clock_t> tps_passe{ };
    clock_t t0 = clock();

    AlgoRSA2D algo2D_Sphe_1(ToreLength_2D, desiredRPhi, exclusionDistance, seed,
        1, "Sphere");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;
    AlgoRSA2D algo2D_Sphe_2(ToreLength_2D, desiredRPhi, exclusionDistance, seed,
        1, "Sphere");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;

    AlgoRSA2D algo2D_Cube_1(ToreLength_2D, desiredRPhi, exclusionDistance, seed,
        1, "Cube");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;
    AlgoRSA2D algo2D_Cube_2(ToreLength_2D, desiredRPhi, exclusionDistance, seed,
        2, "Cube");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;

    AlgoRSA3D algo3D_Cyl_1(ToreLength_3D, desiredRPhi, exclusionDistance, seed,
        1, "Cylinder");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;
    AlgoRSA3D algo3D_Cyl_2(ToreLength_3D, desiredRPhi, exclusionDistance, seed,
        2, "Cylinder");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;

    AlgoRSA3D algo3D_Sph_1(ToreLength_3D, desiredRPhi, exclusionDistance, seed,
        1, "Sphere");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;
    AlgoRSA3D algo3D_Sph_2(ToreLength_3D, desiredRPhi, exclusionDistance, seed,
        2, "Sphere");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;

    Point<3> Tore3D_std({ 64, 64, 64 });
    vector<array<double, 2>> desiredRPhi_std({ array<double, 2> { 1, 1 } });

    AlgoRSA3D algo3D_Tore_1(Tore3D_std, desiredRPhi_std, exclusionDistance,
        seed, 1, "Tore");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;
    AlgoRSA3D algo3D_Tore_2(Tore3D_std, desiredRPhi_std, exclusionDistance,
        seed, 2, "Tore");
    t = clock();
    tps_passe.push_back(t - t0);
    t0 = t;

    for (size_t i = 0; i < tps_passe.size();) {
        cout << "algo 1 : "
            << (static_cast<double>(tps_passe[i])) / CLOCKS_PER_SEC << "s"
            << endl;
        i++;
        cout << "algo 2 : "
            << (static_cast<double>(tps_passe[i])) / CLOCKS_PER_SEC << "s"
            << endl;
        i++;
    }

    //
    cout << "algo           " << "Elapsed CPU time: " << " | nb particles"
        << endl << " 2D, 1st algo :"
        << (static_cast<double>(tps_passe[0])) / CLOCKS_PER_SEC << "s |"
        << algo2D_Sphe_1.getSpheres().size() << endl
        << " 2D, 2nd algo :"
        << (static_cast<double>(tps_passe[1])) / CLOCKS_PER_SEC << "s | "
        << algo2D_Sphe_2.getSpheres().size() << endl;

    cout << " 3D, 1st algo : "
        << (static_cast<double>(tps_passe[8])) / CLOCKS_PER_SEC << "s | "
        << algo3D_Sph_1.getSpheres().size() << endl
        << " 3D, 2nd algo : "
        << (static_cast<double>(tps_passe[9])) / CLOCKS_PER_SEC << "s | "
        << algo3D_Sph_2.getSpheres().size() << endl;

    //
    cout << "Fin algo" << endl;

    /*	cout << "---------------------------------------------------" << endl;
     cout << "algo3D_1" << endl;
     cout << "\n " << algo3D_1.verifySphere() << endl;
     cout << "---------------------------------------------------" << endl;
     cout << "algo3D_2" << endl;
     cout << "\n " << algo3D_2.verifySphere() << endl;
     cout << "---------------------------------------------------" << endl;
     cout << "algo2D_1" << endl;
     cout << "\n " << algo2D_1.verifySphere() << endl;
     cout << "---------------------------------------------------" << endl;
     cout << "algo2D_2" << endl;
     cout << "\n " << algo2D_2.verifySphere() << endl;
     cout << "---------------------------------------------------" << endl;*/

    algo3D_Cyl_2.printDump("sortie_3D_1.dump");
    algo3D_Sph_2.printDump("sortie_3D_2.dump");
    algo2D_Sphe_2.printDump("sortie_2D_1.dump");
    algo2D_Cube_2.printDump("sortie_2D_2.dump");

    algo3D_Cyl_2.printCSV("sortie_csv.csv");
    /*
    algo 1 : 0.014s
    algo 2 : 0.014s
    algo 1 : 0.019s
    algo 2 : 0.01s
    algo 1 : 20.868s
    algo 2 : 13.692s
    algo 1 : 15.066s
    algo 2 : 8.965s
    algo 1 : 2.81s
    algo 2 : 2.867s
     */

    cerr << "END OF VARIOUS TESTS" << endl;
    cerr << "---------------" << endl;
}

inline void TestWP(string nameShape) {
    clock_t t_WP = clock();
    AlgoWP<3> algo{ };
    algo.setExclusionDistance(0.0);
    algo.setBigShape(vector<double> { 10., 10., 10. }, nameShape);
    vector<array<double, 2>> desiredRPhi_2 = { array<double, 2> { 1., 0.50 } };
    algo.setRadiusGenerator(desiredRPhi_2);

    algo.proceed(2, 1);
    clock_t tps_passe_WP = clock() - t_WP;

    cout << "----------------" << endl;
    cout << "Test AlgoWP for " << nameShape << endl;
    cout << (static_cast<double>(tps_passe_WP)) / CLOCKS_PER_SEC << "s |"
        << endl;
    cout << algo.verifySphere() << endl;
    cout << algo.getSpheres().size() << endl;
    algo.printDump(nameShape + "WP_algo.dump");
    cout << "----------------" << endl;
}

inline void TestBool(string nameShape) {
    clock_t t_WP = clock();
    AlgoBool<3> algo{ };
    algo.setExclusionDistance(0.0);
    algo.setBigShape(vector<double> { 10., 10., 10. }, nameShape);
    vector<array<double, 2>> desiredRPhi_2 = { array<double, 2> { 1., 0.50 } };
    algo.setRadiusGenerator(desiredRPhi_2);

    algo.proceed(2, 1);
    clock_t tps_passe_WP = clock() - t_WP;

    cout << "----------------" << endl;
    cout << "Test AlgoWP for " << nameShape << endl;
    cout << (static_cast<double>(tps_passe_WP)) / CLOCKS_PER_SEC << "s |"
        << endl;
    cout << algo.verifySphere() << endl;
    cout << algo.getSpheres().size() << endl;
    algo.printDump(nameShape + "Bool_algo.dump");
    cout << "----------------" << endl;
}

inline void MultiTestBool() {
    cerr << "---------------" << endl;
    cerr << "BEGINNING OF VARIOUS TESTS for Bool" << endl;
    TestBool("Tore");
    TestBool("Cube");
    TestBool("Sphere");
    TestBool("Cylinder");
    cerr << "END OF VARIOUS TESTS for Bool" << endl;
    cerr << "---------------" << endl;
}

inline void MultiTestWP() {
    cerr << "---------------" << endl;
    cerr << "BEGINNING OF VARIOUS TESTS for WP" << endl;
    TestWP("Tore");
    TestWP("Cube");
    TestWP("Sphere");
    TestWP("Cylinder");
    cerr << "END OF VARIOUS TESTS for WP" << endl;
    cerr << "---------------" << endl;
    /*
     FOR NUM_METHOD_OVERLAP = 1;
     ---------------
     BEGINNING OF VARIOUS TESTS for WP
     ----------------------------
     Stopped because no more spheres should be placed.
     Fraction volumique : 0.199989
     Not fully packed
     Nb of spheres : 7639
     ----------------------------
     WPGrid<DIM>::removeOverlap : Number of iterations = 16
     WPGrid<DIM>::removeOverlap : Number of iterations = 41
     WPGrid<DIM>::removeOverlap : Number of iterations = 89
     ----------------
     Test AlgoWP for Tore
     1.701s |
     1
     7639
     ----------------
     ----------------------------
     Stopped because no more spheres should be placed.
     Fraction volumique : 0.199989
     Not fully packed
     Nb of spheres : 7639
     ----------------------------
     WPGrid<DIM>::removeOverlap : Number of iterations = 19
     WPGrid<DIM>::removeOverlap : Number of iterations = 44
     WPGrid<DIM>::removeOverlap : Number of iterations = 136
     ----------------
     Test AlgoWP for Cube
     1.792s |
     1
     7639
     ----------------
     ----------------------------
     Stopped because no more spheres should be placed.
     Fraction volumique : 0.2
     Not fully packed
     Nb of spheres : 4000
     ----------------------------
     WPGrid<DIM>::removeOverlap : Number of iterations = 19
     WPGrid<DIM>::removeOverlap : Number of iterations = 47
     WPGrid<DIM>::removeOverlap : Number of iterations = 119
     ----------------
     Test AlgoWP for Sphere
     1.005s |
     1
     4000
     ----------------
     ----------------------------
     Stopped because no more spheres should be placed.
     Fraction volumique : 0.2
     Not fully packed
     Nb of spheres : 6000
     ----------------------------
     WPGrid<DIM>::removeOverlap : Number of iterations = 21
     WPGrid<DIM>::removeOverlap : Number of iterations = 41
     WPGrid<DIM>::removeOverlap : Number of iterations = 151
     ----------------
     Test AlgoWP for Cylinder
     1.39s |
     1
     6000
     ----------------
     END OF VARIOUS TESTS for WP
     ---------------

     ||||||||||||||||||||||||||||||||||||||||||||||||||||||||
     FOR NUM_METHOD_OVERLAP = 2;
     Test on the 2021/05/10
     ---------------
     BEGINNING OF VARIOUS TESTS for WP
     ----------------------------
     Stopped because no more spheres should be placed.
     Fraction volumique : 0.199989
     Not fully packed
     Nb of spheres : 7639
     ----------------------------
     WPGrid<DIM>::removeOverlap : Number of iterations = 20
     WPGrid<DIM>::removeOverlap : Number of iterations = 41
     WPGrid<DIM>::removeOverlap : Number of iterations = 93
     ----------------
     Test AlgoWP for Tore
     1.758s |
     1
     7639
     ----------------
     ----------------------------
     Stopped because no more spheres should be placed.
     Fraction volumique : 0.199989
     Not fully packed
     Nb of spheres : 7639
     ----------------------------
     WPGrid<DIM>::removeOverlap : Number of iterations = 21
     WPGrid<DIM>::removeOverlap : Number of iterations = 52
     WPGrid<DIM>::removeOverlap : Number of iterations = 145
     ----------------
     Test AlgoWP for Cube
     1.964s |
     1
     7639
     ----------------
     ----------------------------
     Stopped because no more spheres should be placed.
     Fraction volumique : 0.2
     Not fully packed
     Nb of spheres : 4000
     ----------------------------
     WPGrid<DIM>::removeOverlap : Number of iterations = 20
     WPGrid<DIM>::removeOverlap : Number of iterations = 55
     WPGrid<DIM>::removeOverlap : Number of iterations = 134
     ----------------
     Test AlgoWP for Sphere
     0.975s |
     1
     4000
     ----------------
     ----------------------------
     Stopped because no more spheres should be placed.
     Fraction volumique : 0.2
     Not fully packed
     Nb of spheres : 6000
     ----------------------------
     WPGrid<DIM>::removeOverlap : Number of iterations = 22
     WPGrid<DIM>::removeOverlap : Number of iterations = 49
     WPGrid<DIM>::removeOverlap : Number of iterations = 130
     ----------------
     Test AlgoWP for Cylinder
     1.677s |
     1
     6000
     ----------------
     END OF VARIOUS TESTS for WP
     ---------------
     */
}

inline void volInter() {
    // see https://mathworld.wolfram.com/Sphere-SphereIntersection.html
    // see https://mathworld.wolfram.com/Circle-CircleIntersection.html
    Point<3> p1{ 0, 0, 0 };
    Point<3> p2{ 0, 0, 0.69459 };
    AmbiantSpace::Tore<3> tore(Point<3> { 10, 10, 10 });
    double res = sphereTools::volInter(Sphere<3>(p1, 1, 0), Sphere<3>(p2, 1, 0),
        &tore);
    if (abs(res / (m_PI * 4. / 3.) - 0.5) > 0.001) {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
    Point<2> p1_{ 0, 0 };
    Point<2> p2_{ 0, 0.8079455 };
    AmbiantSpace::Tore<2> tore_(Point<2> { 10, 10 });
    double rayon = 1;
    res = sphereTools::volInter(Sphere<2>(p1_, rayon, 0),
        Sphere<2>(p2_, rayon, 0), &tore_);
    if (abs(res / m_PI - 0.5) > 0.001) {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

inline void Interface() {
    algoSpheres::TypeAlgo typeAlgo = algoSpheres::TypeAlgo::RSA;
    constexpr unsigned short DIM = 3;
    Point<DIM> L = { 10., 20., 30. };
    auto nameShape = AmbiantSpace::NameShape::Sphere;
    unsigned seed = 12;
    vector<array<double, 2>> desiredRPhi = { array<double, 2> { 1, 0.2 }, array<
            double, 2> { 2., 0.1 } };
    vector<PhaseType> tabPhases = { 0, 1 };
    double mindist = 0.05;
    algoSpheres::throwSpheres<DIM>(typeAlgo, nameShape, L, seed, desiredRPhi,
        tabPhases, mindist);
}

inline void NamePhases() {
    string nameShape = "Tore";
    AlgoWP<3> algo{ };
    algo.setExclusionDistance(0.0);
    algo.setBigShape(vector<double> { 25., 25., 25. }, nameShape);
    vector<array<double, 2>> desiredRPhi_2 = { array<double, 2> { 5., 0.15 },
            array<double, 2> { 4., 0.05 } };
    vector<PhaseType> tabPhases = { 5, 36 };
    map<PhaseType, string> dico = { { 5, "ph1" }, { 36, "ph2" } };
    algo.setNamePhase(dico);
    algo.setRadiusGenerator(desiredRPhi_2, tabPhases);
    algo.proceed(2, 1);
    algo.printPos(nameShape + "_lb_test.pos");
}

inline void ExclusionDistance() {
    AlgoWP<3> algo{ };
    double L = 100;
    double bdryEDist = 8;
    double eDist = 4;
    double radius = 5;
    algo.setExclusionDistance(eDist);
    algo.setBoundaryExclusionDistance(bdryEDist);
    //
    algo.setBigShape(vector<double> { L, L, L }, "Cube");
    vector<array<double, 2>> desiredRPhi_2 = { array<double, 2> { radius,
            0.13404129 } };
    algo.setRadiusGenerator(desiredRPhi_2);
    //
    algo.proceed(2, 1);
    algo.printPos("WP_algo.pos");

    array<array<double, 2>, 3> posMax;
    for (size_t i = 0; i < 3; i++) {
        posMax[i][0] = numeric_limits<double>::max();
        posMax[i][1] = 0;
    }

    for (auto sph : algo.getSpheres()) {
        for (size_t i = 0; i < 3; i++) {
            posMax[i][0] = min(posMax[i][0], sph.center[i]);
            posMax[i][1] = max(posMax[i][1], sph.center[i]);
        }
    }

    for (size_t j = 0; j < 2; j++) {
        for (size_t i = 0; i < 3; i++) {
            cerr << posMax[i][j] << endl;
        }
    }
    for (size_t i = 0; i < 3; i++) {
        if (posMax[i][0] < bdryEDist + radius
            or posMax[i][1] > L - (bdryEDist + radius)) {
            throw runtime_error("Boundary distance not respected!");
        }
    }
}

inline void PickOnSphere() {
    size_t seed = 1;
    mt19937 engine(seed);
    size_t i_max = 200;
    double L_0 = 4.;
    constexpr unsigned short DIM = 3;
    Point<DIM> L = create_array<DIM>(L_0);
    vector<Sphere<DIM>> theSpheres{};
    for (size_t i = 0; i < i_max; i++) {
        Sphere<DIM> sph{};
        sph.center = randomShooter::pickOnSphere<DIM>(engine);
        for (auto& coord : sph.center) {
            coord += 0.5 * L_0;
        }
        sph.radius = 0.05;
        theSpheres.push_back(sph);
    }
    SphereManipulator<DIM> sphManip(theSpheres, L);
    sphManip.printDump("UniformDistribution_on_sphere.dump");
}

template<unsigned short DIM>
inline void FillMaxRSA() {
    cerr << "Testing " << __PRETTY_FUNCTION__ << endl;
    Point<DIM> L{};
    if constexpr (DIM == 3) {
        L = Point<DIM>{ 10., 5., 9. };
    } else if constexpr (DIM == 2) {
        L = Point<DIM>{ 10., 5. };
    }
    size_t NbSpheres = 2000;
    unsigned seed = 423;
    double mindist = 0.05;
    auto nameShape = AmbiantSpace::NameShape::Tore;
    algoSpheres::fillMaxRSA<DIM>(nameShape, L, NbSpheres, seed, mindist);
}

inline void FillMaxRSA_2() {
    cerr << "Begin Testing " << __PRETTY_FUNCTION__ << endl;
    auto sv = 0;  // Random number generations seed for seperating spheres positionning
    auto N = 20;  // Nombre de grains
    double r = 150;  // Rayons des grains
    double l2D = r * sqrt(m_PI * N / 0.54689);  // Dimension du VER 2D (carrÃ©)
    Point<2> L = { l2D, l2D };
    double mindist = 0;
    auto theSpheres = algoSpheres::fillMaxRSA<2>(AmbiantSpace::NameShape::Tore, L, N, sv, mindist);
    cerr << "End Testing " << __PRETTY_FUNCTION__ << endl;
}

inline void TestPerf_singleTest(int num_threads, double L_0) {
    constexpr unsigned short DIM = 3;
    Point<DIM> L = { L_0, L_0, L_0 };
    omp_set_num_threads(num_threads);
    auto spheres = algoSpheres::throwSpheres<DIM>(algoSpheres::TypeAlgo::WP, AmbiantSpace::NameShape::Tore, L, 0, vector<array<double, 2>>{ { 1, 0.5 }}, { 1 }, 0.);
}

inline void TestPerf() {
    vector<double> L_0_list = { 10, 20, 30, 40, 50, 75, 100 };
    vector<int> nbThreads = { 1, 2, 4, 8, 16, 32 };
    for (auto L_0 : L_0_list) {
        for (auto num_threads : nbThreads) {
            TestPerf_singleTest(num_threads, L_0);
        }
    }
}

inline void GeomTests() {
    {
        // triangle area
        Point<3> x0 = { 0, 0, 0 };
        Point<3> x1 = { 3, 0, 0 };
        Point<3> x2 = { 0, 4, 0 };
        if (abs(geomTools::area::triangle<3>(x0, x1, x2) - 6) > 0.000000001) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Triangle or edge lengths 3-4-5 of incorrect area!");
        }
    }
    {
        // ellipse area
        double radius = 2.4;
        Ellipse<2> ellipse(Sphere<2>({ 0, 0 }, radius, 0));
        for (int i = 0; i < 4; i++) {
            double theta = i * 0.5 * M_PI;
            if (abs(ellipseAux::computeAngularArea(ellipse, theta)
                - 0.5 * theta * radius * radius) > 0.000001 * radius * radius) {
                cerr << i << endl;
                cerr << theta << endl;
                cerr << ellipseAux::computeAngularArea(ellipse, theta) << endl;
                cerr << 0.5 * theta * radius * radius << endl;
                cerr << __PRETTY_FUNCTION__ << endl;
                throw runtime_error("Parts of sphere of incorrect area!");
            }
        }
        Point<2> x1 = radius * Point<2>{ cos(0.25 * M_PI), sin(0.25 * M_PI) };
        Point<2> x2 = radius * Point<2> { cos(0.25 * M_PI), -sin(0.25 * M_PI) };
        if (abs(ellipseAux::computeChordArea(ellipse, x1, x2)
            - radius * radius * (0.75 * M_PI + 0.5)) > 0.000001 * radius * radius) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Chord of incorrect area #1!");
        }
        if (abs(ellipseAux::computeChordArea(ellipse, x2, x1)
            - radius * radius * (0.25 * M_PI - 0.5)) > 0.000001 * radius * radius) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Chord of incorrect area #2!");
        }
        ellipse = Ellipse<2>({ 1, 1 }, { Point<2>{1, 0.}, Point<2>{0, 2.} });
        if (abs(ellipseAux::computeAngularArea(ellipse, M_PI) - 0.5 * ellipse.volume()) > 0.000001) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Incorrect computation of the ellipse area");
        }
    }
}

}  // namespace algoRSA_aux_test
}  // namespace sac_de_billes


