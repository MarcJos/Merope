//! Copyright : see license.txt
//!
//! \brief


#include "Mesh/TestMesh.hxx"
#include "../../AlgoPacking/include/Interface.hxx"
#include "Mesh/GmshWriter.hxx"
#include "MultiInclusions/MultiInclusions.hxx"
#include "MicroToMesh/MeshGenerator.hxx"
#include "Voronoi/VoroToMeshGraph.hxx"


namespace merope {

void testMesh::test1() {
    using namespace mesh;
    using namespace mesh::meshStructure;
    using namespace mesh::geoObjects;

    constexpr short DIM = 3;

    Point<DIM>              L{ 1, 1, 1 };
    // 1
    GeoPoint<DIM> pt1(4, Point<DIM>{0, 0, 0});
    GeoPoint<DIM> pt2(1, Point<DIM>{1, 0, 0});
    GeoPoint<DIM> pt3(2, Point<DIM>{0, 1, 0});
    GeoPoint<DIM> pt4(3, Point<DIM>{0, 0, 1});
    GeoPoint<DIM> pt5(100, Point<DIM>{0, 0, 1.0000001});
    vector<GeoPoint<DIM>>  vGeo{ pt1, pt2, pt3, pt4, pt5 };
    // 2
    Edge e1(4, { 4, 1 });
    Edge e2(5, { 4, 2 });
    Edge e3(6, { 4, 3 });
    Edge e4(7, { 1, 2 });
    Edge e5(8, { 1, 3 });
    Edge e6(9, { 2, 3 });
    Edge e7(101, { 100, 1 });
    vector<Edge>            vecEdge{ e1, e2, e3, e4, e5, e6, e7 };
    // 3
    CurveLoop c1(10, { 5, -7, -4 });
    CurveLoop c2(11, { 4, 8, -6 });
    CurveLoop c3(12, { -5, 6, -9 });
    CurveLoop c4(13, { 7, 9, -8 });
    vector<CurveLoop>       vecCurveLoop{ c1, c2, c3, c4 };
    // 4
    Surface s1(1, { 10 });
    Surface s2(2, { 11 });
    Surface s3(3, { 12 });
    Surface s4(4, { 13 });
    // 5
    vector<Surface>         vecSurface{ s1, s2, s3, s4 };
    // 6
    SurfaceLoop sL(1, { 1, 2, 3, 4 });
    vector<SurfaceLoop>     vecSurfaceLoop{ sL };
    // 7
    Solid so(1, { 1 });
    vector<Solid>           vecSolid{ so };
    //

    mesh::meshStructure::VoroMesh_UnStructureData<DIM> rawData{ L, vGeo, vecEdge, vecCurveLoop, vecSurface, vecSurfaceLoop, vecSolid, {}, {} };
    //!
    mesh::meshStructure::VoroMesh_Periodic<DIM> geoPerStruct(rawData, false);
    //!
    string nameFile = "OutputFile.geo";
    ofstream f{ nameFile };
    mesh::gmsh_writer::write<mesh::gmsh_writer::MeshMethod::native_gmsh>(geoPerStruct, f);
}

void testMesh::test2() {
    constexpr unsigned short DIM = 3;
    double l = 1.;
    double mindist = 0.01 * l;
    Point<DIM> L = { l, l, l };
    auto typeAlgo{ algoSpheres::TypeAlgo::RSA };
    auto nameShape = AmbiantSpace::NameShape::Tore;
    vector<array<double, 2>> desiredRPhi = { {0.2 * l, 1} };
    vector<PhaseType> tabPhases = { };
    long seed = 0;
    vector<Sphere<DIM>> theSpheres = algoSpheres::throwSpheres<3>(typeAlgo, nameShape, L, seed, desiredRPhi, tabPhases, mindist);


    // DEBUG TEST
    // NE PAS EFFACER : PB VORO++
    ofstream myFile("spheres.txt");
    for (size_t i = 0; i < theSpheres.size(); i++) {
        const auto& sphere = theSpheres[i];
        myFile << i + 1 << " " << sphere.center[0] << " " << sphere.center[1] << " " << sphere.center[2] << " " << sphere.radius << endl;
    }
    voro::container_poly con(0, L[0], 0, L[1], 0, L[2], 6, 6, 6, true, true, true, 8);

    // Import the monodisperse test packing and output the Voronoi
    // tessellation in gnuplot and POV-Ray formats.
    con.import("spheres.txt");
    // Do a custom output routine to store the number of vertices, edges,
    // and faces of each Voronoi cell
    con.print_custom("ID=%i, pos=(%x,%y,%z), vertices=%w, edges=%g, faces=%s", "packing.custom1");
    con.print_custom("ID=%i, vertices = %p", "packing.custom2");
    // NE PAS EFFACER : PB VORO++

    // DEBUG TEST

    microToMesh::voroTranslater::VoroToMeshGraph voroTranslater(L, theSpheres);
    mesh::meshStructure::VoroMesh_UnStructureData<3> rawMeshData = voroTranslater.getMeshData();
    mesh::meshStructure::VoroMesh_Periodic<DIM> geoPerStructure(rawMeshData, true);
    ofstream f1("Output_0.geo");
    //geoPerStructure.print(cerr);
    mesh::gmsh_writer::write<mesh::gmsh_writer::MeshMethod::native_gmsh>(geoPerStructure, f1);

    ofstream f2("Output_0_enveloppe.geo");
    geoPerStructure.restrictEnveloppe();
    mesh::gmsh_writer::write<mesh::gmsh_writer::MeshMethod::native_gmsh>(geoPerStructure, f2);
}

void testMesh::test3() {
    constexpr unsigned short DIM = 3;
    array<double, DIM> L = create_array<DIM>(1.);
    size_t nbSpheres = 4;
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, nbSpheres, 0, 0.05);
    for (auto& sphere : theSpheres) {
        sphere.phase = 2;
    }
    //////////////////
    SphereInclusions<DIM> sphInc{};
    sphInc.setLength(L);
    sphInc.setSpheres(theSpheres);

    MultiInclusions<DIM> mi{};
    mi.setInclusions(sphInc);
    mi.setMatrixPhase(1);
    // mi.addLayer(mi.getAllIdentifiers(), 2, 0.1);

    mesh::generator::MeshGenerator meshGenerator{};
    meshGenerator.setMeshOrder(2);
    meshGenerator.setMeshSize(0.05);
    meshGenerator.setMultiInclusions(mi);
    meshGenerator.write("Spheres_mesh.geo");
}

void testMesh::test4() {
    cerr << __PRETTY_FUNCTION__ << endl;
    constexpr unsigned short DIM = 3;
    array<double, DIM> L = create_array<DIM>(1.);
    size_t nbSpheres = 5;
    auto theSpheres = algoSpheres::fillMaxRSA<DIM>(AmbiantSpace::NameShape::Tore, L, nbSpheres, 0, 0.05);
    for (auto& sphere : theSpheres) {
        sphere.phase = 2;
    }
    //////////////////
    LaguerreTess<DIM> sphInc(L, theSpheres);

    MultiInclusions<DIM> mi{};
    mi.setInclusions(sphInc);
    mi.addLayer(mi.getAllIdentifiers(), 3, 0.05);

    mesh::generator::MeshGenerator meshGenerator{};
    meshGenerator.setMeshOrder(2);
    meshGenerator.setMeshSize(0.02);
    meshGenerator.setMultiInclusions(mi);
    meshGenerator.write("LaguerreTess_mesh.geo");
}

void testMesh::test5() {
    Point<3> L = { 2, 2, 2 };
    unsigned seed = 0;
    vector<array<double, 2>> desiredRPhi = { {0.5, 1} };
    vector<PhaseType> tabPhases = { 1 };
    auto spheres = sac_de_billes::algoSpheres::throwSpheres<3>(sac_de_billes::algoSpheres::TypeAlgo::RSA, sac_de_billes::AmbiantSpace::NameShape::Tore, L, seed, desiredRPhi, tabPhases, 0.);
    // print
    SphereInclusions<3> mIncl{};
    mIncl.setLength(L);
    mIncl.setSpheres(spheres);
    mIncl.printVER("Seeds.txt");
    // print
    LaguerreTess<3> polyCrystal(L, spheres);
    MultiInclusions<3> multiInclusions{};
    multiInclusions.setInclusions(polyCrystal);
    auto new_phases = multiInclusions.getAllIdentifiers();
    for (auto& ph : new_phases) ph += 2;
    multiInclusions.setMatrixPhase(1);
    multiInclusions.changePhase(multiInclusions.getAllIdentifiers(), new_phases);


    mesh::generator::MeshGenerator meshGenerator{};
    meshGenerator.setMeshOrder(1);
    meshGenerator.setMeshSize(0.025);
    meshGenerator.setMultiInclusions(multiInclusions);
    meshGenerator.setAdimMergeDistance(1.e-5);
    meshGenerator.write("mmm_mesh.geo");
}

void testMesh::test6(int nb_sph) {
    Point<3> L = { 1, 1, 1 };
    int nbSpheres = nb_sph;
    double distMin = 0.05;
    int randomSeed = 0;
    auto theSpheres = sac_de_billes::algoSpheres::fillMaxRSA<3>(AmbiantSpace::NameShape::Tore, L, nbSpheres, randomSeed, distMin);

    for (auto& sphere : theSpheres) {
        sphere.phase = 2;
    }

    merope::SphereInclusions<3>  sphInc{};
    sphInc.setLength(L);
    sphInc.setSpheres(theSpheres);

    MultiInclusions<3> mi{};
    mi.setInclusions(sphInc);
    mi.setMatrixPhase(1);
    mi.addLayer(mi.getAllIdentifiers(), 3, 0.05);
    mi.addLayer(mi.getAllIdentifiers(), 4, 0.05);


    mesh::generator::MeshGenerator meshGenerator{};
    meshGenerator.setMeshOrder(2);
    meshGenerator.setMeshSize(0.025);
    meshGenerator.setMultiInclusions(mi);
    meshGenerator.do_not_mesh({ 2 });
    meshGenerator.set_nameOutput({ "spheres.vtk" });
    meshGenerator.write("sphere_mesh.geo");
}

}  // namespace merope

