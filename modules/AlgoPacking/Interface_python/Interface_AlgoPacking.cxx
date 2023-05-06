//! \file Interface.cxx
//! \date 26 avr. 2021
//!
//! \brief 
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include<iostream>
#include<omp.h>

#include "../src/AlgoRSA.hxx"
#include "../src/AlgoPacking.cxx"
#include "../src/SphereManipulator.hxx"
#include "../src/AlgoWP.hxx"
#include "../src/AlgoBool.hxx"
#include "../src/Interface.hxx"
#include "../src/AlgoNames.hxx"
#include "../src/GlobalShape.hxx"

namespace py = pybind11;

using namespace sac_de_billes;


int get_num_threads(){
int result = 0;
#pragma omp parallel
#pragma omp master
    result = omp_get_num_threads();
return result;
}

void createModule(py::module_& rsa){

py::class_<Sphere<3>>(rsa,"Sphere_3D")
    .def(py::init<>())
    .def(py::init<const Point<3> &, double, int>())
    .def_readwrite("center", &Sphere<3>::center)
    .def_readwrite("radius", &Sphere<3>::radius)
    .def_readwrite("phase", &Sphere<3>::phase)
    ;

py::class_<Sphere<2>>(rsa,"Sphere_2D")
    .def(py::init<>())
    .def(py::init<const Point<2> &, double, int>())
    .def_readwrite("center", &Sphere<2>::center)
    .def_readwrite("radius", &Sphere<2>::radius)
    .def_readwrite("phase", &Sphere<2>::phase)
    ;

py::class_<AlgoRSA3D>(rsa,"AlgoRSA3D") // for backward compatibility
	.def(py::init<>())
	.def(py::init<std::array<double, 3>, std::vector<std::array<double, 2> >, double, unsigned, unsigned short, std::string>())
	.def(py::init<std::array<double, 3>, std::vector<std::array<double, 2> >, double, unsigned, unsigned short, AmbiantSpace::NameShape>())
	;

py::class_<AlgoRSA<3>>(rsa,"AlgoRSA_3D")
	.def(py::init<>())
	.def("setNamePhase", &AlgoRSA<3>::setNamePhase)
	.def("proceed", &AlgoRSA<3>::proceed)
	.def("printDump", &AlgoRSA<3>::printDump)
	.def("printCSV", &AlgoRSA<3>::printCSV)
	.def("printPos", &AlgoRSA<3>::printPos)
	.def("getPlacedSpheres", &AlgoRSA<3>::getPlacedSpheres_nonT)
	.def("getSpheres", &AlgoRSA<3>::getSpheres)
	.def("getPhases", &AlgoRSA<3>::getPhases)
	.def("getLength", &AlgoRSA<3>::getLength)
	.def("verifySphere", &AlgoRSA<3>::verifySphere)
	.def("setBigShape", static_cast<void (AlgoRSA<3>::*)(vector<double>, std::string)>(&AlgoRSA<3>::setBigShape))
	.def("setBigShape", static_cast<void (AlgoRSA<3>::*)(vector<double>, AmbiantSpace::NameShape)>(&AlgoRSA<3>::setBigShape))
	.def("setExclusionDistance", &AlgoRSA<3>::setExclusionDistance)
	.def("setBoundaryExclusionDistance", &AlgoRSA<3>::setBoundaryExclusionDistance)
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<3>::*)(std::vector<std::array<double, 2>>, std::vector<PhaseType>)>(&AlgoRSA<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<3>::*)(std::vector<double>, std::vector<PhaseType>)>(&AlgoRSA<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<3>::*)(std::vector<std::array<double, 2>>)>(&AlgoRSA<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<3>::*)(std::vector<double>)>(&AlgoRSA<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<3>::*)(std::string nameFile)>(&AlgoRSA<3>::setRadiusGenerator))
	;

py::class_<AlgoRSA<2>>(rsa,"AlgoRSA_2D")
	.def(py::init<>())
	.def("setNamePhase", &AlgoRSA<2>::setNamePhase)
	.def("proceed", &AlgoRSA<2>::proceed)
	.def("printDump", &AlgoRSA<2>::printDump)
	.def("printCSV", &AlgoRSA<2>::printCSV)
	.def("getPlacedSpheres", &AlgoRSA<2>::getPlacedSpheres_nonT)
	.def("getSpheres", &AlgoRSA<2>::getSpheres)
	.def("getPhases", &AlgoRSA<2>::getPhases)
	.def("getLength", &AlgoRSA<2>::getLength)
	.def("verifySphere", &AlgoRSA<2>::verifySphere)
	.def("setBigShape", static_cast<void (AlgoRSA<2>::*)(vector<double>, std::string)> (&AlgoRSA<2>::setBigShape))
	.def("setBigShape", static_cast<void (AlgoRSA<2>::*)(vector<double>, AmbiantSpace::NameShape)> (&AlgoRSA<2>::setBigShape))
	.def("setExclusionDistance", &AlgoRSA<2>::setExclusionDistance)
	.def("setBoundaryExclusionDistance", &AlgoRSA<2>::setBoundaryExclusionDistance)
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<2>::*)(std::vector<std::array<double, 2>>, std::vector<PhaseType>)>(&AlgoRSA<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<2>::*)(std::vector<double>, std::vector<PhaseType>)>(&AlgoRSA<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<2>::*)(std::vector<std::array<double, 2>>)>(&AlgoRSA<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<2>::*)(std::vector<double>)>(&AlgoRSA<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoRSA<2>::*)(std::string nameFile)>(&AlgoRSA<2>::setRadiusGenerator))
	;

py::class_<AlgoBool<3>>(rsa,"AlgoBool_3D")
	.def(py::init<>())
	.def("setNamePhase", &AlgoBool<3>::setNamePhase)
	.def("proceed", &AlgoBool<3>::proceed)
	.def("printDump", &AlgoBool<3>::printDump)
	.def("printCSV", &AlgoBool<3>::printCSV)
	.def("printPos", &AlgoBool<3>::printPos)
	.def("getPlacedSpheres", &AlgoBool<3>::getPlacedSpheres_nonT)
	.def("getSpheres", &AlgoBool<3>::getSpheres)
	.def("getPhases", &AlgoBool<3>::getPhases)
	.def("getLength", &AlgoBool<3>::getLength)
	.def("verifySphere", &AlgoBool<3>::verifySphere)
	.def("setBigShape", static_cast<void (AlgoBool<3>::*)(vector<double>, std::string)>(&AlgoBool<3>::setBigShape))
	.def("setBigShape", static_cast<void (AlgoBool<3>::*)(vector<double>, AmbiantSpace::NameShape)>(&AlgoBool<3>::setBigShape))
	.def("setExclusionDistance", &AlgoBool<3>::setExclusionDistance)
	.def("setBoundaryExclusionDistance", &AlgoBool<3>::setBoundaryExclusionDistance)
	.def("setRadiusGenerator", static_cast<void (AlgoBool<3>::*)(std::vector<std::array<double, 2>>, std::vector<PhaseType>)>(&AlgoBool<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoBool<3>::*)(std::vector<double>, std::vector<PhaseType>)>(&AlgoBool<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoBool<3>::*)(std::vector<std::array<double, 2>>)>(&AlgoBool<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoBool<3>::*)(std::vector<double>)>(&AlgoBool<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoBool<3>::*)(std::string nameFile)>(&AlgoBool<3>::setRadiusGenerator))
	;

py::class_<AlgoBool<2>>(rsa,"AlgoBool_2D")
	.def(py::init<>())
	.def("setNamePhase", &AlgoBool<2>::setNamePhase)
	.def("proceed", &AlgoBool<2>::proceed)
	.def("printDump", &AlgoBool<2>::printDump)
	.def("printCSV", &AlgoBool<2>::printCSV)
	.def("getPlacedSpheres", &AlgoBool<2>::getPlacedSpheres_nonT)
	.def("getSpheres", &AlgoBool<2>::getSpheres)
	.def("getPhases", &AlgoBool<2>::getPhases)
	.def("getLength", &AlgoBool<2>::getLength)
	.def("verifySphere", &AlgoBool<2>::verifySphere)
	.def("setBigShape", static_cast<void (AlgoBool<2>::*)(vector<double>, std::string)> (&AlgoBool<2>::setBigShape))
	.def("setBigShape", static_cast<void (AlgoBool<2>::*)(vector<double>, AmbiantSpace::NameShape)> (&AlgoBool<2>::setBigShape))
	.def("setExclusionDistance", &AlgoBool<2>::setExclusionDistance)
	.def("setBoundaryExclusionDistance", &AlgoBool<2>::setBoundaryExclusionDistance)
	.def("setRadiusGenerator", static_cast<void (AlgoBool<2>::*)(std::vector<std::array<double, 2>>, std::vector<PhaseType>)>(&AlgoBool<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoBool<2>::*)(std::vector<double>, std::vector<PhaseType>)>(&AlgoBool<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoBool<2>::*)(std::vector<std::array<double, 2>>)>(&AlgoBool<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoBool<2>::*)(std::vector<double>)>(&AlgoBool<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoBool<2>::*)(std::string nameFile)>(&AlgoBool<2>::setRadiusGenerator))
	;

py::class_<AlgoWP<3>>(rsa,"AlgoWP_3D")
	.def(py::init<>())
	.def("setNamePhase", &AlgoWP<3>::setNamePhase)
	.def("proceed", &AlgoWP<3>::proceed)
	.def("printDump", &AlgoWP<3>::printDump)
	.def("printCSV", &AlgoWP<3>::printCSV)
	.def("printPos", &AlgoWP<3>::printPos)
	.def("getPlacedSpheres", &AlgoWP<3>::getPlacedSpheres_nonT)
	.def("getSpheres", &AlgoWP<3>::getSpheres)
	.def("getPhases", &AlgoWP<3>::getPhases)
	.def("getLength", &AlgoWP<3>::getLength)
	.def("verifySphere", &AlgoWP<3>::verifySphere)
	.def("setBigShape", static_cast<void (AlgoWP<3>::*)(vector<double>, std::string)>(&AlgoWP<3>::setBigShape))
	.def("setBigShape", static_cast<void (AlgoWP<3>::*)(vector<double>, AmbiantSpace::NameShape)>(&AlgoWP<3>::setBigShape))
	.def("setExclusionDistance", &AlgoWP<3>::setExclusionDistance)
	.def("setBoundaryExclusionDistance", &AlgoWP<3>::setBoundaryExclusionDistance)
	.def("setRadiusGenerator", static_cast<void (AlgoWP<3>::*)(std::vector<std::array<double, 2>>, std::vector<PhaseType>)>(&AlgoWP<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoWP<3>::*)(std::vector<double>, std::vector<PhaseType>)>(&AlgoWP<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoWP<3>::*)(std::vector<std::array<double, 2>>)>(&AlgoWP<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoWP<3>::*)(std::vector<double>)>(&AlgoWP<3>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoWP<3>::*)(std::string nameFile)>(&AlgoWP<3>::setRadiusGenerator))
	;

py::class_<AlgoWP<2>>(rsa,"AlgoWP_2D")
	.def(py::init<>())
	.def("setNamePhase", &AlgoWP<2>::setNamePhase)
	.def("proceed", &AlgoWP<2>::proceed)
	.def("printDump", &AlgoWP<2>::printDump)
	.def("printCSV", &AlgoWP<2>::printCSV)
	.def("getPlacedSpheres", &AlgoWP<2>::getPlacedSpheres_nonT)
	.def("getSpheres", &AlgoWP<2>::getSpheres)
	.def("getPhases", &AlgoWP<2>::getPhases)
	.def("getLength", &AlgoWP<2>::getLength)
	.def("verifySphere", &AlgoWP<2>::verifySphere)
	.def("setBigShape", static_cast<void (AlgoWP<2>::*)(vector<double>, std::string)> (&AlgoWP<2>::setBigShape))
	.def("setBigShape", static_cast<void (AlgoWP<2>::*)(vector<double>, AmbiantSpace::NameShape)> (&AlgoWP<2>::setBigShape))
	.def("setExclusionDistance", &AlgoWP<2>::setExclusionDistance)
	.def("setBoundaryExclusionDistance", &AlgoWP<2>::setBoundaryExclusionDistance)
	.def("setRadiusGenerator", static_cast<void (AlgoWP<2>::*)(std::vector<std::array<double, 2>>, std::vector<PhaseType>)>(&AlgoWP<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoWP<2>::*)(std::vector<double>, std::vector<PhaseType>)>(&AlgoWP<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoWP<2>::*)(std::vector<std::array<double, 2>>)>(&AlgoWP<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoWP<2>::*)(std::vector<double>)>(&AlgoWP<2>::setRadiusGenerator))
	.def("setRadiusGenerator", static_cast<void (AlgoWP<2>::*)(std::string nameFile)>(&AlgoWP<2>::setRadiusGenerator))
	;

py::class_<SphereManipulator<2>>(rsa,"SphereManipulator_2D")
	.def("setNamePhase", &SphereManipulator<2>::setNamePhase)
	.def(py::init<const AlgoInterface<2,algoWP_aux::AlgoWP_Template<2>> &>())
	.def(py::init<const AlgoInterface<2,algoRSA_aux::AlgoRSA_Template<2>> &>())
	.def(py::init<std::vector<SimpleSphereFormat>, std::array<double,2>>())
	.def(py::init<std::vector<Sphere<2>>, std::array<double, 2>>())
	.def("printVER", &SphereManipulator<2>::printVER)
	.def("printDump", &SphereManipulator<2>::printDump)
	.def("printCSV", &SphereManipulator<2>::printCSV)
	.def("getPlacedSpheres", &SphereManipulator<2>::getPlacedSpheres_nonT)
	.def("getSpheres", &SphereManipulator<2>::getSpheres)
	.def("getLength", &SphereManipulator<2>::getLength)
	.def("translate", &SphereManipulator<2>::translate)
	.def("translate_back", &SphereManipulator<2>::translate_back)
	.def("distance_cube_spheres", &SphereManipulator<2>::distance_cube_spheres)
	.def("random_search", &SphereManipulator<2>::random_search)
	.def("upper_bound_on_best_distmin", &SphereManipulator<2>::upper_bound_on_best_distmin)
	;

py::class_<SphereManipulator<3>>(rsa,"SphereManipulator_3D")
	.def("setNamePhase", &SphereManipulator<3>::setNamePhase)
	.def(py::init<const AlgoInterface<3,algoWP_aux::AlgoWP_Template<3>> &>())
	.def(py::init<const AlgoInterface<3,algoRSA_aux::AlgoRSA_Template<3>> &>())
	.def(py::init<std::vector<SimpleSphereFormat>, std::array<double,3>>())
	.def(py::init<std::vector<Sphere<3>>, std::array<double, 3>>())
	.def("printVER", &SphereManipulator<3>::printVER)
	.def("printDump", &SphereManipulator<3>::printDump)
	.def("printCSV", &SphereManipulator<3>::printCSV)
	.def("printPos", &SphereManipulator<3>::printPos)
	.def("getPlacedSpheres", &SphereManipulator<3>::getPlacedSpheres_nonT)
	.def("getSpheres", &SphereManipulator<3>::getSpheres)
	.def("getLength", &SphereManipulator<3>::getLength)
	.def("translate", &SphereManipulator<3>::translate)
	.def("translate_back", &SphereManipulator<3>::translate_back)
	.def("distance_cube_spheres", &SphereManipulator<3>::distance_cube_spheres)
	.def("random_search", &SphereManipulator<3>::random_search)
	.def("upper_bound_on_best_distmin", &SphereManipulator<3>::upper_bound_on_best_distmin)
	;


// Interface.hxx


rsa.def("fillMaxRSA_3D", &algoSpheres::fillMaxRSA<3>);
rsa.def("fillMaxRSA_2D", &algoSpheres::fillMaxRSA<2>);

//openmp
rsa.def("omp_set_num_threads", &omp_set_num_threads);
rsa.def("omp_get_num_threads", &get_num_threads);

//rsa.def("fromCumHisto_3D", &algoSpheres::fromCumHisto<3>);
//rsa.def("fromCumHisto_2D", &algoSpheres::fromCumHisto<2>);

rsa.def("throwSpheres_3D", static_cast<std::vector<Sphere<3>> (*)(const algoSpheres::TypeAlgo &,
        AmbiantSpace::NameShape, std::array<double, 3>, unsigned, std::vector<std::array<double, 2>>, vector<PhaseType>, double)>(&algoSpheres::throwSpheres<3>));
rsa.def("throwSpheres_2D", static_cast<std::vector<Sphere<2>> (*)(const algoSpheres::TypeAlgo &,
        AmbiantSpace::NameShape, std::array<double, 2>, unsigned, std::vector<std::array<double, 2>>, vector<PhaseType>, double)>(&algoSpheres::throwSpheres<2>));
rsa.def("throwSpheres_3D", static_cast<std::vector<Sphere<3>> (*)(const algoSpheres::TypeAlgo &,
        AmbiantSpace::NameShape, std::array<double, 3>, unsigned, std::vector<double>, vector<PhaseType>, double)>(&algoSpheres::throwSpheres<3>));
rsa.def("throwSpheres_2D", static_cast<std::vector<Sphere<2>> (*)(const algoSpheres::TypeAlgo &,
        AmbiantSpace::NameShape, std::array<double, 2>, unsigned, std::vector<double>, vector<PhaseType>, double)>(&algoSpheres::throwSpheres<2>));

// enum
py::enum_<algoSpheres::TypeAlgo>(rsa,"TypeAlgo")
    .value("RSA", algoSpheres::TypeAlgo::RSA)
    .value("BOOL", algoSpheres::TypeAlgo::BOOL)
    .value("WP", algoSpheres::TypeAlgo::WP)
    .export_values()
    ;

py::enum_<AmbiantSpace::NameShape>(rsa,"NameShape")
    .value("Tore", AmbiantSpace::NameShape::Tore)
    .value("Cube", AmbiantSpace::NameShape::Cube)
    .value("Sphere", AmbiantSpace::NameShape::Sphere)
    .value("Cylinder", AmbiantSpace::NameShape::Cylinder)
    .export_values()
    ;
}



PYBIND11_MODULE(py_rsa, rsa){
	createModule(rsa);
}

PYBIND11_MODULE(sac_de_billes, sac_de_billes_){
	createModule(sac_de_billes_);
}
