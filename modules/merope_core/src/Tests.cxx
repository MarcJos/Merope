//! Copyright : see license.txt
//!
//! \brief

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Test/Tests.hxx"
#include "../../AlgoPacking/src/GlobalShape.hxx"
#include "../../AlgoPacking/src/Interface.hxx"
#include "../../AlgoPacking/src/For_testing.hxx"
#include "MultiInclusions/MultiInclusions.hxx"
#include "MultiInclusions/SphereInclusions.hxx"
#include "Obsolete_MesoStructure/InterfaceStructure.hxx"
#include "Obsolete_MesoStructure/MicroType.hxx"
#include "Voxellation/Voxellation.hxx"
#include "VTKinout/VTKStream.hxx"
#include "Parallelism/SetNbOfThreads.hxx"
#include "AlgoLaguerre/Optimize_LaguerreTess.hxx"
#include "AlgoLaguerre/MakeCentroidal.hxx"

#include "MeropeNamespace.hxx"


namespace merope {

void Tests::polyCrystal11() {
    auto sv = 0; // Random number generations seed for seperating spheres positionning
    auto N = 20; // Nb spheres
    double l3D = 10; // RVE dimensions
    constexpr unsigned short DIM = 3;
    array<double, DIM> L = { l3D, l3D, l3D };
    double mindist = 0;
    // Discretization
    size_t n3D = 64;
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, N, sv, mindist);
    double e = 0.125; // width of the layer

    LaguerreTess<DIM> lag(L, theSpheres);
    lag.setAspRatio(array<double, DIM> { 2, 1., 1.});
    MultiInclusions<DIM> mi{};
    mi.setInclusions(lag);

    mi.addLayer(mi.getAllIdentifiers(), 2, e);


    vector<array<double, 2>> desiredRPhi = { {0.4,  0.1} };
    auto spherePores = algoSpheres::throwSpheres<DIM>(algoSpheres::TypeAlgo::BOOL,
        AmbiantSpace::NameShape::Tore, L, sv, desiredRPhi,
        vector<PhaseType> { 1 });
    cerr << spherePores.size() << endl;
    MultiInclusions<DIM> mi2{};
    SphereInclusions<DIM> si{};
    si.setLength(L);
    si.setSpheres(spherePores);
    mi2.setInclusions(si);
    //

    Structure<DIM> structure(mi, mi2, std::map<PhaseType, PhaseType>{ { 2, 3 }});

    for (auto i : structure.getAllPhases()) {
        cerr << "Phase =" << i << endl;
    }

    vox::Voxellation<DIM> voxGrid{ structure };

    voxGrid.proceed(array<size_t, DIM> { n3D, n3D, n3D });

    voxGrid.printFile("poro3D_1.vtk", "Coeffs.txt");

    //
    //
    Grid g{};
}

void Tests::polyCrystal0() {
    auto sv = 0; // Random number generations seed for seperating spheres positionning
    auto N = 20; // Nb spheres
    double l3D = 10; // RVE dimensions
    constexpr unsigned short DIM = 3;
    array<double, DIM> L = { l3D, l3D, l3D };
    double mindist = 0;
    // Discretization
    size_t n3D = 64;
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, N, sv, mindist);

    double e = 0.125; // width of the layer
    // Voronoi avec porositÃ©s

    InterfaceStructure<3> v3D{};
    v3D.setLength(L);
    v3D.mainInclusions.setTypeCrystal(TypeCrystal::Laguerre);
    v3D.mainInclusions.setAspRatio(array<double, 3> { 2, 1., 1.});
    v3D.mainInclusions.setSpheres(theSpheres);

    v3D.setErosionWidth(e);

    vector<array<double, 2>> desiredRPhi = { {0.4,  0.1} };
    auto spherePores = algoSpheres::throwSpheres<3>(algoSpheres::TypeAlgo::BOOL, AmbiantSpace::NameShape::Tore, L, sv, desiredRPhi, vector<PhaseType>{1});
    v3D.secdInclusions.setSpheres(spherePores);
    v3D.setColorization(ColorMaterialID::Erode2Mat);

    auto structure_0 = v3D.build();
    // Grille
    vox::Voxellation<3> voxGrid_0{ structure_0 };
    voxGrid_0.setPureCoeffs(vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21});
    voxGrid_0.proceed(array<size_t, 3> {n3D, n3D, n3D});
    voxGrid_0.printFile("poro3D_2_0.vtk", "Coeffs_3D_0.txt");

    v3D.setColorization(ColorMaterialID::Erode3Mat);
    auto structure = v3D.build();
    // Grille
    vox::Voxellation<3> voxGrid{ structure };
    voxGrid.setPureCoeffs(vector<double>{1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 15, 16, 17, 18, 19, 20, 21});
    voxGrid.proceed(array<size_t, 3> {n3D, n3D, n3D});
    voxGrid.printFile("poro3D_2.vtk", "Coeffs_3D.txt");

    InterfaceStructure<3> vNew{};
    vNew.setLength(L);
    vNew.mainInclusions.setTypeCrystal(TypeCrystal::Spheres);
    vNew.mainInclusions.setSpheres(spherePores);
    vNew.setColorization(ColorMaterialID::Erode2Mat);
    vNew.setInnerPhase(1);
    vox::Voxellation<3> newVoxGrid{ vNew.build() };
    newVoxGrid.proceed(array<size_t, 3> {n3D, n3D, n3D});
    newVoxGrid.printFile("poro3D_spheres_2.vtk", "Coeffs_3D.txt");


    vox::symmetrize<3>("poro3D_spheres_2.vtk", "poro3D_spheres_2_sym.vtk", array<size_t, 3>{0, 1, 2});
}



void Tests::polyCrystal10() {
    auto sv = 0; // Random number generations seed for seperating spheres positionning
    auto N = 150; // Nb spheres
    double l3D = 10; // RVE dimensions
    constexpr unsigned short DIM = 2;
    array<double, DIM> L = { l3D, l3D };
    double mindist = 0;
    // Discretization
    size_t n2D = 64;
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, N, sv, mindist);

    double e = 0.125; // width of the layer

    LaguerreTess<DIM> lag(L, theSpheres);
    lag.setAspRatio(array<double, 2> { 2, 1. });
    MultiInclusions<DIM> mi{};
    mi.setInclusions(lag);

    mi.addLayer(mi.getAllIdentifiers(), 2, e);

    vector<array<double, 2>> desiredRPhi = { { 1, 0.1 } };
    auto spherePores = algoSpheres::throwSpheres<2>(algoSpheres::TypeAlgo::BOOL,
        AmbiantSpace::NameShape::Tore, L, sv, desiredRPhi,
        vector<PhaseType> { 1 });
    cerr << spherePores.size() << endl;
    MultiInclusions<2> mi2{};
    SphereInclusions<DIM> si{};
    si.setLength(L);
    si.setSpheres(spherePores);
    mi2.setInclusions(si);

    Structure<DIM> structure(mi, mi2, std::map<PhaseType, PhaseType>{ { 2, 3 }});

    vox::Voxellation<2> voxGrid{ structure };
    voxGrid.setPureCoeffs(vector<double>{1, 2, 3, 4});

    voxGrid.proceed(array<size_t, 2> { n2D, n2D });

    voxGrid.printFile("poro2D_1.vtk", "Coeffs.txt");
}


void Tests::polyCrystal1() {
    auto sv = 0; // Random number generations seed for seperating spheres positionning
    auto N = 150; // Nb spheres
    double l3D = 10; // RVE dimensions
    constexpr unsigned short DIM = 2;
    array<double, DIM> L = { l3D, l3D };
    double mindist = 0;
    // Discretization
    size_t n2D = 256;
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, N, sv, mindist);

    double e = 0.125; // width of the layer

    ofstream file("fractionsBOOL.txt");
    // Voronoi avec porositÃ©s

    InterfaceStructure<2> v2D{};
    v2D.setLength(L);
    v2D.mainInclusions.setTypeCrystal(TypeCrystal::Laguerre);
    v2D.mainInclusions.setAspRatio(array<double, 2> { 2, 1.});
    v2D.mainInclusions.setSpheres(theSpheres);

    v2D.setErosionWidth(e);

    vector<array<double, 2>> desiredRPhi = { {0.1, 0.01} };
    auto spherePores = algoSpheres::throwSpheres<2>(algoSpheres::TypeAlgo::BOOL, AmbiantSpace::NameShape::Tore, L, sv, desiredRPhi, vector<PhaseType>{1});
    v2D.secdInclusions.setSpheres(spherePores);
    v2D.setColorization(ColorMaterialID::Erode3Mat);

    auto structure = v2D.build();
    // Grille
    vox::Voxellation<2> voxGrid{ structure };

    voxGrid.proceed(array<size_t, 2> {n2D, n2D});

    voxGrid.printFile("poro2D_2.vtk", "Coeffs.txt");
}


void Tests::polyCrystal9() {
    auto sv = 0; // Random number generations seed for seperating spheres positionning
    auto N = 150; // Nb spheres
    double l3D = 10; // RVE dimensions
    constexpr unsigned short DIM = 3;
    array<double, DIM> L = { l3D, l3D, l3D };
    double mindist = 0;
    // Discretization
    size_t n3D = 256;
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, N, sv, mindist);

    double e = 0.125; // width of the layer

    InterfaceStructure<DIM> v3D{};
    v3D.setLength(L);
    v3D.mainInclusions.setTypeCrystal(TypeCrystal::Laguerre);
    v3D.mainInclusions.setAspRatio(array<double, DIM> { 2, 1., 3.});
    v3D.mainInclusions.setSpheres(theSpheres);

    v3D.setErosionWidth(e);

    vector<array<double, 2>> desiredRPhi = { {0.25, 0.1} };
    auto spherePores = algoSpheres::throwSpheres<DIM>(algoSpheres::TypeAlgo::BOOL, AmbiantSpace::NameShape::Tore, L, sv, desiredRPhi, vector<PhaseType>{1});
    v3D.secdInclusions.setSpheres(spherePores);
    v3D.setColorization(ColorMaterialID::Erode3Mat);

    auto structure = v3D.build();
    // Grille
    vox::Voxellation<DIM> voxGrid{ structure };

    voxGrid.proceed(array<size_t, DIM> {n3D, n3D, n3D});

    voxGrid.printFile("poro3D-gdVolume.vtk", "Coeffs.txt");
}

template<unsigned short DIM>
inline void auxi_polyCrystal_fit_volumes(string fileVTK) {
    auto sv = 0; // Random number generations seed for seperating spheres positionning
    auto N = 200; // Nb spheres
    double l3D = 10; // RVE dimensions
    array<double, DIM> L = create_array<DIM>(l3D);
    double mindist = 0;
    // Discretization
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, N, sv, mindist);
    for (size_t i = 0; i < theSpheres.size(); i++) {
        theSpheres[i].phase = i;
    }
    /////////////////////////////////////
    vector<double> volumeFractions(theSpheres.size(), auxi_function::productOf<double>(L) / theSpheres.size());
    optimizeLaguerreTess::algo_fit_volumes<DIM> algo(L, theSpheres, volumeFractions);
    clock_t t0 = clock();
    /////////////////////////////
    ///// Important lines
    /////////////////////////////
    algo.proceed(1e-5, 300, true);
    auto vols = algo.getCurrentVolumes();
    for (size_t i = 0; i < vols.size(); i++) {
        std::cerr << vols[i] << " , ";
    }
    std::cerr << endl;
    auto new_center_tessels = algo.getCenterTessels();
    std::cerr << "Nb of spheres" << vols.size() << " = " << new_center_tessels.size() << endl;
    /////////////////////////////
    /////////////////////////////
    double total_time = (static_cast<double>(clock() - t0)) / CLOCKS_PER_SEC;

    cerr << "############" << endl;
    cerr << "############" << endl;
    std::cerr << "Here is the total time :" << total_time << endl;
    cerr << "############" << endl;
    cerr << "############" << endl;

    /////////////////////////////// PRINT
    LaguerreTess<DIM> lag(L, new_center_tessels);
    MultiInclusions<DIM> mi{};
    mi.setInclusions(lag);
    //////////////////////////////////////
    vox::Voxellation<DIM> voxellation{ mi };
    array<size_t, DIM> nbVox = create_array<DIM, size_t>(64);
    voxellation.proceed(nbVox);
    //////////////////////////////////////
    string fileCoeff = "Coeffs.txt";
    voxellation.printFile(fileVTK, fileCoeff);
}

void Tests::polyCrystal_fit_volumes_3D(string fileVTK) {
    auxi_polyCrystal_fit_volumes<3>(fileVTK);
}

void Tests::polyCrystal_fit_volumes_2D(string fileVTK) {
    auxi_polyCrystal_fit_volumes<2>(fileVTK);
}

void Tests::polyCrystalCentroidal(bool use_acceleration) {
    auto sv = 0; // Random number generations seed for seperating spheres positionning
    auto N = 4000; // Nb spheres
    double l3D = 10; // RVE dimensions
    constexpr unsigned short DIM = 3;
    array<double, DIM> L = { l3D, l3D, l3D };
    double mindist = 0;
    // Discretization
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, N, sv, mindist);
    optimizeLaguerreTess::makeCentroidal<DIM, 4>(L, theSpheres, 100, 1e-5, use_acceleration, true);
    std::cerr << theSpheres.size() << endl;
}


void Tests::polyCrystal8() {
    auto sv = 0; // Random number generations seed for seperating spheres positionning
    auto N = 150; // Nb spheres
    double l3D = 10; // RVE dimensions
    constexpr unsigned short DIM = 3;
    array<double, DIM> L = { l3D, l3D, l3D };
    double mindist = 0;
    // Discretization
    size_t n3D = 256;
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, N, sv, mindist);

    double e = 0.125; // width of the layer

    InterfaceStructure<DIM> v3D{};
    v3D.setLength(L);
    v3D.mainInclusions.setTypeCrystal(TypeCrystal::Laguerre);
    v3D.mainInclusions.setAspRatio(array<double, DIM> { 2, 1., 3.});
    v3D.mainInclusions.setSpheres(theSpheres);

    v3D.setErosionWidth(e);
    v3D.setColorization(ColorMaterialID::Poly);

    auto structure = v3D.build();
    // Grid
    vox::Voxellation<DIM> voxGrid{ structure };

    voxGrid.proceed(array<size_t, DIM> {n3D, n3D, n3D});

    voxGrid.printFile("poly3D-gdVolume.vtk", "Coeffs.txt");
}

void Tests::spheres1() {
    constexpr unsigned short DIM = 3;
    array<double, DIM> L = create_array<DIM>(10.);
    array<size_t, DIM> nbVox = create_array<DIM, size_t>(64);
    size_t nbSpheres = 1500;
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, nbSpheres, 0, 0);
    //////////////////
    SphereInclusions<DIM> sphInc{};
    sphInc.setLength(L);
    sphInc.setSpheres(theSpheres);

    MultiInclusions<DIM> mi{};
    mi.setInclusions(sphInc);
    mi.addLayer(mi.getAllIdentifiers(), 1, 0.1);

    vox::Voxellation<DIM> voxellation{ mi };
    voxellation.proceed(nbVox);
    string fileVTK = "Spheres.vtk";
    string fileCoeff = "Coeffs.txt";
    voxellation.printFile(fileVTK, fileCoeff);
}

void Tests::extraction() {
    constexpr unsigned short DIM = 3;
    Point<DIM> L = { 10, 10, 10 };
    size_t n3D = 256;
    array<size_t, DIM> nVox = { n3D, n3D, n3D };
    array<size_t, DIM> nVoxMin = { 100, 100, 100 };
    array<size_t, DIM> nVoxMax = { 101, 150, 200 };

    SphereInclusions<DIM> sphIncl2;
    sphIncl2.setLength(L);
    sphIncl2.fromHisto(0, algoSpheres::TypeAlgo::RSA, 0., vector<array<double, 2>>{ { 0.25, 1 }}, vector<long>{1});

    MultiInclusions<DIM> multiInclusions2{};
    multiInclusions2.setInclusions(sphIncl2);

    vox::Voxellation<DIM> grid1(multiInclusions2);
    grid1.setPureCoeffs(vector<double>{0, 1});
    grid1.setVoxelRule(vox::VoxelRule::Center);
    grid1.setHomogRule(homogenization::Rule::Largest);
    grid1.proceed(nVox, nVoxMin, nVoxMax);
    grid1.printFile("Zone_Inclusions_extraction.vtk", "Coeffs.txt");

    SphereInclusions<DIM> sphIncl{};
    sphIncl.setLength(L);
    sphIncl.fromHisto(0, algoSpheres::TypeAlgo::RSA, 0., vector<array<double, 2>>{ { 2., 1. }}, vector<long>{1});
    sphIncl.printVER("Seeds.txt");

    LaguerreTess<DIM> polyCrystal(L, sphIncl.getSpheres());
    MultiInclusions<DIM> multiInclusions{};
    multiInclusions.setInclusions(polyCrystal);
    size_t N = multiInclusions.getAllIdentifiers().size();
    multiInclusions.addLayer(multiInclusions.getAllIdentifiers(), N, 0.3);
    vector<long> newId(multiInclusions.getAllIdentifiers().size(), 1);
    multiInclusions.changePhase(multiInclusions.getAllIdentifiers(), newId);

    vector<double> coeffs{};
    for (size_t i = 0; i < N + 1; i++) {
        coeffs.push_back(i);
    }

    vox::Voxellation<DIM> grid2(multiInclusions);
    grid2.setPureCoeffs(coeffs);
    grid2.setVoxelRule(vox::VoxelRule::Center);
    grid2.setHomogRule(homogenization::Rule::Largest);
    grid2.proceed(nVox, nVoxMin, nVoxMax);
    grid2.printFile("Zone_Crystal_extraction.vtk", "Coeffs.txt");

    map<PhaseType, PhaseType> dictionnaire{};
    dictionnaire[N] = N + 1;
    Structure<DIM> structure(multiInclusions, multiInclusions2, dictionnaire);

    vox::Voxellation<DIM> grid3(structure);
    coeffs.push_back(N + 1);
    grid3.setPureCoeffs(coeffs);
    grid3.setVoxelRule(vox::VoxelRule::Center);
    grid3.setHomogRule(homogenization::Rule::Largest);
    grid3.proceed(nVox, nVoxMin, nVoxMax);
    grid3.printFile("Zone_extraction.vtk", "Coeffs.txt");
}

void Tests::testPerf0() {
    const long NB_XP = 50;
    constexpr short DIM = 2;
    array<double, DIM> L = { 20., 20. };
    array<size_t, DIM> nbVox = { 609, 1176 };
    double total_time = 0;
    for (auto i = 0; i < NB_XP; i++) {
        SphereInclusions<DIM> sphIncl{};
        sphIncl.setLength(L);
        sphIncl.fromHisto(i, algoSpheres::TypeAlgo::RSA, 0, vector<array<double, 2>>{ { 0.3, 0.2 }, { 0.3, 0.2 }}, vector<PhaseType>{1, 1});
        ///
        MultiInclusions<DIM> multiInclusions{};
        multiInclusions.setInclusions(sphIncl);
        vox::Voxellation<DIM> grid{ multiInclusions };
        grid.setHomogRule(homogenization::Rule::Voigt);
        grid.setVoxelRule(vox::VoxelRule::Average);
        clock_t t0 = clock();
        grid.proceed(nbVox);
        total_time += (static_cast<double>(clock() - t0)) / CLOCKS_PER_SEC;
        cerr << "Total time = " << total_time << endl;
    }
    for (size_t i = 0; i < 10; i++) {
        cerr << "############" << endl;
    }
    clock_t t0 = clock();
    vector<double> phases(609 * 1176);
    for (size_t i = 0; i < NB_XP; i++) {
        phases[i] = 9.09;
    }
    cerr << "Total time = " << (static_cast<double>(clock() - t0)) / CLOCKS_PER_SEC << endl;
}

void Tests::testFields() {
    constexpr unsigned short DIM = 3;
    Point<DIM> L{ 10, 10, 10 };
    std::function<double(Point<DIM>)> covariance = [](Point<DIM> x) {return exp(-geomTools::normeCarre<DIM>(x));};
    std::function<double(double)> nonlinearFunction = [](double a) {return exp(a);};
    CartesianField<DIM> cartesianField(gaussianField::SimpleGaussianField<DIM>(covariance, nonlinearFunction), L);
    FieldStructure<DIM> fieldStructure(cartesianField);
    vox::Voxellation<DIM> myVox(fieldStructure);
}

void Tests::testFields2() {
    constexpr unsigned short DIM = 3;
    double r_0 = 2;
    double r_1 = 3;
    double r_2 = 4;
    double l_silver = 0.3;
    double phi_1 = 0.2;
    double  phi_2 = 0.3;
    double lambda_0clay = 0.6;
    double lambda_lead = 35.3;
    double lambda_silver = 429;

    array<double, DIM> L = { 50, 50, 50 };
    array<size_t, DIM> N3D = { 128, 128, 128 };
    int seed_sph = 0;
    int seed_Gauss = 1;

    // define gaussian field
    auto covariance = [r_0](Point<DIM> x) {return exp((-x[0] * x[0] - x[1] * x[1] - x[2] * x[2]) / (2 * r_0 * r_0));};
    auto nonLin = [lambda_0clay](double g) {return lambda_0clay * (2 + tanh(g));};
    merope::gaussianField::SimpleGaussianField<DIM> gaussianne(covariance, nonLin);
    gaussianne.seed = seed_Gauss;
    // get spheres
    vector<array<double, 2>> desiredRPhi = { {r_1, phi_1}, {r_2, phi_2} };
    vector<long> tabPhases = { 1, 1 };
    double minDist = 0.;
    auto the_spheres = sac_de_billes::algoSpheres::throwSpheres<3>(sac_de_billes::algoSpheres::TypeAlgo::WP, sac_de_billes::AmbiantSpace::NameShape::Tore, L, seed_sph, desiredRPhi, tabPhases, minDist);
    // spherical inclusions
    merope::SphereInclusions<3> sphIncl{};
    sphIncl.setLength(L);
    sphIncl.setSpheres(the_spheres);
    merope::MultiInclusions<3> mIncl{};
    mIncl.setInclusions(sphIncl);
    auto list_id = mIncl.getAllIdentifiers();
    mIncl.addLayer(list_id, 2, l_silver);
    mIncl.setMatrixPhase(2);
    // mask
    merope::SphereInclusions<3> maskSphIncl{};
    maskSphIncl.setLength(L);
    maskSphIncl.setSpheres(the_spheres);
    merope::MultiInclusions<3> mask{};
    mask.setInclusions(sphIncl);
    mask.setMatrixPhase(0);
    // printStruc
    auto printStruc = [](const auto& struc, const auto& N3D_, const string& name) {
        merope::vox::Voxellation<DIM> voxelTot(struc);
        voxelTot.proceed(N3D_);
        voxelTot.printFile(name + ".vtk", name + ".txt");
        };
    // print Gaussian Field
    merope::CartesianField<DIM> cField_Gauss(gaussianne, L);
    merope::FieldStructure<DIM> struc_gauss(cField_Gauss);
    merope::vox::Voxellation<DIM> vox(struc_gauss);
    vox.proceed(N3D);
    vox.printFile("gauss.vtk", "gauss.txt");
    printStruc(struc_gauss, N3D, "struc_gauss");
    // print mask
    merope::Structure<DIM> structureMask(mask);
    merope::vox::Voxellation<DIM> voxellationMask(structureMask);
    voxellationMask.setVoxelRule(merope::vox::VoxelRule::Average);
    voxellationMask.setHomogRule(merope::homogenization::Rule::Voigt);
    voxellationMask.proceed(N3D);
    voxellationMask.printFile("mask.vtk", "mask.txt");
    auto fieldMask = voxellationMask.getField();
    merope::CartesianField<DIM> cField_mask(fieldMask, L);
    merope::FieldStructure<DIM> struc_mask(cField_mask);
    printStruc(struc_mask, N3D, "struc_mask");
    // print inclusions
    merope::Structure<DIM> structureIncl(mIncl);
    merope::vox::Voxellation<DIM> voxellationIncl(structureIncl);
    voxellationIncl.setVoxelRule(merope::vox::VoxelRule::Average);
    voxellationIncl.setHomogRule(merope::homogenization::Rule::Reuss);
    voxellationIncl.setPureCoeffs({ lambda_0clay, lambda_lead, lambda_silver });
    voxellationIncl.proceed(N3D);
    voxellationIncl.printFile("inclusions.vtk", "inclusions.txt");
    auto fieldInclusions = voxellationIncl.getField();
    merope::CartesianField<DIM> cField_inclusions(fieldInclusions, L);
    merope::FieldStructure<DIM> struc_incl(cField_inclusions);
    printStruc(struc_incl, N3D, "struc_incl");
    // resulting microstructure
    merope::FieldStructure<DIM>fieldStructureTot(struc_gauss, struc_incl, struc_mask);
    merope::vox::Voxellation<DIM> voxelTot(fieldStructureTot);
    voxelTot.setVoxelRule(merope::vox::VoxelRule::Average);
    voxelTot.setHomogRule(merope::homogenization::Rule::Reuss);
    voxelTot.proceed(N3D);
    voxelTot.printFile("totalStruct.vtk", "totalStruct.txt");
}


} // namespace merope

