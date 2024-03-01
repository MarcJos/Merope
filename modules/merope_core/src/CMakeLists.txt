include_directories("${PROJECT_SOURCE_DIR}/modules/merope_core/include")
link_directories("${CMAKE_CURRENT_BINARY_DIR}")

find_package(OpenMP)

add_library(Physics Homogenization.cxx)
install(TARGETS Physics DESTINATION lib)

add_library(VoroInterface VoroInterface.cxx MicroInclusion.cxx)
target_link_libraries(VoroInterface PUBLIC vorolib)
target_include_directories(VoroInterface SYSTEM PUBLIC ../../voro-plus-plus)
install(TARGETS VoroInterface DESTINATION lib)

add_library(GeomVER Rotation.cxx Hyper.cxx SphereSeeds.cxx RecurStructure.cxx) 
target_include_directories(GeomVER PUBLIC ${CMAKE_CURRENT_SOURCE_DIR})
target_link_libraries(GeomVER PUBLIC RSAalgolib VoroInterface)
install(TARGETS GeomVER DESTINATION lib)

add_library(OptiLaguerre Optimize_LaguerreTess.cxx)
target_link_libraries(OptiLaguerre PUBLIC VoroInterface)
install(TARGETS OptiLaguerre DESTINATION lib)

add_library(GeneOrientations geneOrientations.cxx)
target_link_libraries(GeneOrientations PUBLIC GeomVER)
install(TARGETS GeneOrientations DESTINATION lib)

add_library(VTKinout VTKStream.cxx VTKRead.cxx TIFFRead.cxx)
install(TARGETS VTKinout DESTINATION lib)

add_library(FFTGrid OpenmpWrapper.cxx FFTField.cxx FFTScalarField.cxx Geostat.cxx CovSum.cxx GaussianField.cxx)
###########
### Add MKL
	###	target_link_libraries(FFTGrid PUBLIC OpenMP::OpenMP_CXX) -> bad interaction with python lambdas in the gaussian functions
	target_link_libraries(FFTGrid PUBLIC -Wl,--start-group ${MKL_ROOT_LIB}/libmkl_intel_ilp64.a ${MKL_ROOT_LIB}/libmkl_intel_thread.a ${MKL_ROOT_LIB}/libmkl_core.a -Wl,--end-group -lomp5 -lpthread -lm -ldl)
	target_compile_options(FFTGrid PUBLIC -DMKL_ILP64  -m64)
	target_include_directories(FFTGrid PRIVATE "${MKL_INCLUDE_DIRS}")
##########
install(TARGETS FFTGrid DESTINATION lib)

add_library(Grid Grid_VER.cxx Grid.cxx CartesianGrid.cxx GridManipulations.cxx SetNbOfThreads.cxx)
target_link_libraries(Grid PUBLIC RSAalgolib VTKinout GeomVER Physics)
target_link_libraries(Grid PUBLIC OpenMP::OpenMP_CXX)
target_link_libraries(Grid PRIVATE FFTGrid)
install(TARGETS Grid DESTINATION lib)

add_library(Mesh GeoObjects.cxx SameThings.cxx GmshWriter.cxx)
target_link_libraries(Mesh PUBLIC RSAalgolib)
install(TARGETS Mesh DESTINATION lib)

add_library(VoroToMesh VoroToMeshGraph.cxx MeshGenerator.cxx)
target_link_libraries(VoroToMesh PUBLIC Mesh VoroInterface)

add_library(tests Tests.cxx TestMesh.cxx)
target_link_libraries(tests PUBLIC RSAalgolib GeomVER VTKinout Grid VoroInterface Mesh VoroToMesh GeneOrientations OptiLaguerre)
install(TARGETS tests DESTINATION lib)


add_executable(Main_pour_test Z_Main_pour_test.cxx)
target_link_libraries(Main_pour_test PUBLIC tests)
install(TARGETS Main_pour_test DESTINATION bin)

add_executable(Laguerre_opti Z_Main_Laguerre_opti.cxx)
target_link_libraries(Laguerre_opti PUBLIC tests)
install(TARGETS Laguerre_opti DESTINATION bin)


add_executable(Performance Z_Performance.cxx)
target_link_libraries(Performance PUBLIC tests)
install(TARGETS Performance DESTINATION bin)


add_executable(Performance_1 Z_Performance_1.cxx)
target_link_libraries(Performance_1 PUBLIC tests)
install(TARGETS Performance_1 DESTINATION bin)

## See Makefile.am in previous versions -> is no C++
## add_executable(ntfa_evp3d_castem ntfa_evp3d_castem.f)
