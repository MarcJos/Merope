//! Copyright : see license.txt
//!
//! \brief 
//

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl_bind.h>
#include <pybind11/functional.h>
#include <pybind11/numpy.h>
#include<iostream>
#include<functional>

#include "../AlgoPacking/src/AlgoNames.hxx"
#include "../include/Parallelism/SetNbOfThreads.hxx"
#include "../include/MultiInclusions/SphereInclusions.hxx"
#include "../include/Voronoi/VoroInterface.hxx"
#include "../include/MultiInclusions/LaguerreTess.hxx"
#include "../include/MesoStructure/Structure.hxx"
#include "../include/Obsolete_MesoStructure/TypeCrystal.hxx"
#include "../include/Obsolete_MesoStructure/InterfaceStructure.hxx"
#include "../include/Field/GaussianField.hxx"
#include "../include/Field/CovSum.hxx"
#include "../include/Grid/Geostat.hxx"
#include "../include/Grid/Grid_VER.hxx"
#include "../include/Voxellation/Voxellation.hxx"
#include "../include/Voxellation/VoxRecurStructure.hxx"
#include "../include/geneOrientations.hxx"
#include "../include/VTKinout/VTKRead.hxx"
#include "../include/VTKinout/TIFFRead.hxx"
#include "../include/VTKinout/VTKStream.hxx"
#include "../include/Grid/GridTypes.hxx"
#include "../include/MicroToMesh/MeshGenerator.hxx"
#include "../include/MicroInclusion/MicroInclusion.hxx"

namespace py = pybind11;

using namespace merope;
using namespace sac_de_billes;

// Voxellation


// stop the GIL of Python in order to allow Python functions to be used in C++
// beware, this is dangerous if Python creates/deletes objects
class Stop_gil{
public:
    Stop_gil():_save{nullptr}{
        _save = PyEval_SaveThread();
    }
    ~Stop_gil(){
        PyEval_RestoreThread(_save);
    }

private:
    PyThreadState *_save;
};

template<class T>
py::array_t<T> create_matrix(vector<T> myVector){
    return py::array_t<T>(py::buffer_info(myVector.data(), sizeof(T), py::format_descriptor<T>::format(),
            1, std::vector<size_t> {myVector.size()}, std::vector<size_t> {sizeof(T)}));
}

template<class VOX>
py::array_t<double> getField_Numpy(VOX& vox){
    vector<double> result{};
    {// stop the GIL of python just for this computation
        Stop_gil stop_gil{};
        result = vox.getField();
    }
    return create_matrix(result);
}

template<class VOX>
py::array_t<vox::VTK_PHASE> getPhaseGrid(VOX& vox){
    vector<vox::VTK_PHASE> result{};
    {// stop the GIL of python just for this computation
        Stop_gil stop_gil{};
        result = vox.getPhaseGrid();
    }
    return create_matrix(result);
}


void createModule_merope(py::module_& merope) {
//VTKinout
py::class_<VTKstream>(merope,"VTKstream")
    .def(py::init<const char*>())
    .def("close",&VTKstream::close)
    ;

py::class_<VTKRead>(merope,"VTKRead")
    .def(py::init<const char*>())
    ;

py::class_<TIFFRead>(merope,"TIFFRead")
    .def(py::init<const char*>())
    ;

py::class_<vox::GridField<3>>(merope, "GridField_3D")
    ;

py::class_<vox::GridField<2>>(merope, "GridField_2D")
    ;

py::class_<vox::Voxellation<3>>(merope,"Voxellation_3D")
	.def(py::init<const MultiInclusions<3>&>())
	.def(py::init<const Structure<3>&>())
	.def(py::init<const FieldStructure<3>&>())
    .def("proceed", static_cast<void(vox::Voxellation<3>::*)(std::array<size_t, 3>, std::array<size_t, 3>, std::array<size_t, 3>)>(&vox::Voxellation<3>::proceed), py::call_guard<Stop_gil>())
    .def("proceed", static_cast<void(vox::Voxellation<3>::*)(std::array<size_t, 3>)>(&vox::Voxellation<3>::proceed), py::call_guard<Stop_gil>())
    .def("computePhaseGrid", &vox::Voxellation<3>::computePhaseGrid, py::call_guard<Stop_gil>())
	.def("printFile", &vox::Voxellation<3>::printFile, py::call_guard<Stop_gil>())
	.def("getField_Numpy", &getField_Numpy<vox::Voxellation<3>>)
	.def("getField", &vox::Voxellation<3>::getField)
	.def("getPhaseGrid", &getPhaseGrid<vox::Voxellation<3>>)
	.def("setPureCoeffs", &vox::Voxellation<3>::setPureCoeffs)
	.def("setHomogRule",&vox::Voxellation<3>::setHomogRule)
	.def("setVoxelRule",&vox::Voxellation<3>::setVoxelRule)
	;

py::class_<vox::Voxellation<2>>(merope,"Voxellation_2D")
    .def(py::init<const MultiInclusions<2>&>())
    .def(py::init<const Structure<2>&>())
    .def(py::init<const FieldStructure<2>&>())
    .def("proceed", static_cast<void(vox::Voxellation<2>::*)(std::array<size_t, 2>, std::array<size_t, 2>, std::array<size_t, 2>)>(&vox::Voxellation<2>::proceed), py::call_guard<Stop_gil>())
    .def("proceed", static_cast<void(vox::Voxellation<2>::*)(std::array<size_t, 2>)>(&vox::Voxellation<2>::proceed), py::call_guard<Stop_gil>())
    .def("printFile", &vox::Voxellation<2>::printFile, py::call_guard<Stop_gil>())
    .def("setPureCoeffs", &vox::Voxellation<2>::setPureCoeffs, py::call_guard<Stop_gil>())
    .def("computePhaseGrid", &vox::Voxellation<2>::computePhaseGrid, py::call_guard<Stop_gil>())
    .def("getField_Numpy", &getField_Numpy<vox::Voxellation<2>>)
    .def("getField", &vox::Voxellation<2>::getField)
    .def("getPhaseGrid", &getPhaseGrid<vox::Voxellation<2>>)
//    .def("extractGrid",&vox::Voxellation<2>::extractGrid)
    .def("setHomogRule",&vox::Voxellation<2>::setHomogRule)
    .def("setVoxelRule",&vox::Voxellation<2>::setVoxelRule)
    ;

merope.def("get_linear_index_2D", &vox::auxi::get_linear_index<2, array<size_t, 2>>);
merope.def("get_linear_index_3D", &vox::auxi::get_linear_index<3, array<size_t, 3>>);

py::class_<Grid_VER>(merope,"Grid")
    .def(py::init<>())
    .def(py::init<VTKRead &>())
    .def("DivNormFactor",&Grid_VER::DivNormFactor)
    .def("toVTKCELL",&Grid_VER::toVTKCELL)
    .def("VTKheaderCELL",&Grid_VER::VTKheaderCELL)
    .def("phaseFraction",&Grid_VER::phaseFraction)
  ;

//! Microstructure
py::class_<SphereInclusions<3>>(merope,"SphereInclusions_3D")
    .def(py::init<>())
    .def("setSpheres", static_cast<void(SphereInclusions<3>::*)(const std::vector<Sphere<3>>)>(&SphereInclusions<3>::setSpheres))
    .def("fromFile",static_cast<void(SphereInclusions<3>::*)(std::string)>(&SphereInclusions<3>::fromFile))
    .def("fromFile",static_cast<void(SphereInclusions<3>::*)(std::string, std::string)>(&SphereInclusions<3>::fromFile))
    .def("setLength",  &SphereInclusions<3>::setLength)
    .def("getLength", &SphereInclusions<3>::getLength)
    .def("printDump",  &SphereInclusions<3>::printDump)
    .def("fromHisto", &SphereInclusions<3>::fromHisto)
    .def("printFracVol",&SphereInclusions<3>::printFracVol)
    .def("printDump",&SphereInclusions<3>::printDump)
    .def("printPos",&SphereInclusions<3>::printPos)
    .def("printCSV",&SphereInclusions<3>::printCSV_space)
    .def("printVER", &SphereInclusions<3>::printVER)
    .def("getSpheres", &SphereInclusions<3>::getSpheres)
    .def("translate", &SphereInclusions<3>::translate)
    ;

py::class_<SphereInclusions<2>>(merope,"SphereInclusions_2D")
    .def(py::init<>())
    .def("setSpheres", static_cast<void(SphereInclusions<2>::*)(const std::vector<Sphere<2>>)>(&SphereInclusions<2>::setSpheres))
    .def("fromFile",static_cast<void(SphereInclusions<2>::*)(std::string)>(&SphereInclusions<2>::fromFile))
    .def("fromFile",static_cast<void(SphereInclusions<2>::*)(std::string, std::string)>(&SphereInclusions<2>::fromFile))
    .def("setLength",  &SphereInclusions<2>::setLength)
    .def("getLength",  &SphereInclusions<2>::getLength)
    .def("printDump",  &SphereInclusions<2>::printDump)
    .def("fromHisto", &SphereInclusions<2>::fromHisto)
    .def("printFracVol",&SphereInclusions<2>::printFracVol)
    .def("printDump",&SphereInclusions<2>::printDump)
    .def("printPos",&SphereInclusions<2>::printPos)
    .def("printCSV",&SphereInclusions<2>::printCSV_space)
    .def("printVER", &SphereInclusions<2>::printVER)
    .def("getSpheres", &SphereInclusions<2>::getSpheres)
    .def("translate", &SphereInclusions<2>::translate)
    ;

py::class_<PolyInclusions<3>>(merope, "PolyInclusions_3D")
    .def(py::init<>())
    .def("setLength",  &PolyInclusions<3>::setLength)
    .def("addInclusion", &PolyInclusions<3>::addInclusion)
    .def("setInclusions", &PolyInclusions<3>::setInclusions)
    ;

py::class_<PolyInclusions<2>>(merope, "PolyInclusions_2D")
    .def(py::init<>())
    .def("setLength",  &PolyInclusions<2>::setLength)
    .def("addInclusion", &PolyInclusions<2>::addInclusion)
    .def("setInclusions", &PolyInclusions<2>::setInclusions)
    ;

py::class_<SpheroPolyhedronInclusions<3>>(merope, "SpheroPolyInclusions_3D")
    .def(py::init<>())
    .def("setLength",  &SpheroPolyhedronInclusions<3>::setLength)
    .def("addInclusion", &SpheroPolyhedronInclusions<3>::addInclusion)
    .def("setInclusions", &SpheroPolyhedronInclusions<3>::setInclusions)
    ;

py::class_<SpheroPolyhedronInclusions<2>>(merope, "SpheroPolyInclusions_2D")
    .def(py::init<>())
    .def("setLength",  &SpheroPolyhedronInclusions<2>::setLength)
    .def("addInclusion", &SpheroPolyhedronInclusions<2>::addInclusion)
    .def("setInclusions", &SpheroPolyhedronInclusions<2>::setInclusions)
    ;


py::class_<MultiInclusions<3>>(merope,"MultiInclusions_3D")
	.def(py::init<>())
	.def("setInclusions", static_cast<void (MultiInclusions<3>::*)(LaguerreTess<3>)>(&MultiInclusions<3>::setInclusions))
	.def("setInclusions", static_cast<void (MultiInclusions<3>::*)(const PolyInclusions<3>&)>(&MultiInclusions<3>::setInclusions))
	.def("setInclusions", static_cast<void (MultiInclusions<3>::*)(const SphereInclusions<3> &)>(&MultiInclusions<3>::setInclusions))
    .def("setInclusions", static_cast<void (MultiInclusions<3>::*)(const EllipseInclusions<3>&)>(&MultiInclusions<3>::setInclusions))
    .def("setInclusions", static_cast<void (MultiInclusions<3>::*)(const std::vector<smallShape::SpheroPolyhedronInc<3>>&, Point<3>)>(&MultiInclusions<3>::setInclusions))
    .def("setInclusions", static_cast<void (MultiInclusions<3>::*)(const SpheroPolyhedronInclusions<3>&)>(&MultiInclusions<3>::setInclusions))
	.def("addLayer", static_cast<void (MultiInclusions<3>::*)(const std::vector<Identifier>& identifiers, PhaseType newPhase, double width)>(&MultiInclusions<3>::addLayer))
	.def("addLayer", static_cast<void (MultiInclusions<3>::*)(const std::vector<Identifier>& identifiers, const std::vector<PhaseType>& newPhase, const std::vector<double>& width)>(&MultiInclusions<3>::addLayer))
	.def("changePhase",&MultiInclusions<3>::changePhase)
	.def("getAllIdentifiers",static_cast<vector<Identifier>(MultiInclusions<3>::*)() const>(&MultiInclusions<3>::getAllIdentifiers))
	.def("getAllPhases",static_cast<std::vector<PhaseType>  (MultiInclusions<3>::*)() const>(&MultiInclusions<3>::getAllPhases))
	.def("getAllPhases",static_cast<std::vector<PhaseType> (MultiInclusions<3>::*)(size_t) const>(&MultiInclusions<3>::getAllPhases))
	.def("getAllCenters",&MultiInclusions<3>::getAllCenters)
	.def("getIdentifiers", &MultiInclusions<3>::getIdentifiers)
	.def("setMatrixPhase", &MultiInclusions<3>::setMatrixPhase)
    //.def("translate", &MultiInclusions<3>::translate)
	//.def("setInnerShapes", static_cast<void (MultiInclusions<3>::*)(const vector<smallShape::ConvexPolyhedron<DIM>>&)>(&MultiInclusions<3>::setInnerShapes))
	//.def("setInnerShapes", static_cast<void (MultiInclusions<3>::*)(const vector<smallShape::SphereInc<DIM>>&)>(&MultiInclusions<3>::setInnerShapes))
	//.def("", &MultiInclusions<3>::)
	;

py::class_<MultiInclusions<2>>(merope,"MultiInclusions_2D")
    .def(py::init<>())
    .def("setInclusions", static_cast<void (MultiInclusions<2>::*)(LaguerreTess<2>)>(&MultiInclusions<2>::setInclusions))
    .def("setInclusions", static_cast<void (MultiInclusions<2>::*)(const PolyInclusions<2>&)>(&MultiInclusions<2>::setInclusions))
    .def("setInclusions", static_cast<void (MultiInclusions<2>::*)(const SphereInclusions<2> &)>(&MultiInclusions<2>::setInclusions))
    .def("setInclusions", static_cast<void (MultiInclusions<2>::*)(const EllipseInclusions<2>&)>(&MultiInclusions<2>::setInclusions))
    .def("setInclusions", static_cast<void (MultiInclusions<2>::*)(const std::vector<smallShape::SpheroPolyhedronInc<2>>&, Point<2>)>(&MultiInclusions<2>::setInclusions))
    .def("addLayer", static_cast<void (MultiInclusions<2>::*)(const std::vector<Identifier>& identifiers, PhaseType newPhase, double width)>(&MultiInclusions<2>::addLayer))
    .def("setInclusions", static_cast<void (MultiInclusions<2>::*)(const SpheroPolyhedronInclusions<2>&)>(&MultiInclusions<2>::setInclusions))
    .def("addLayer", static_cast<void (MultiInclusions<2>::*)(const std::vector<Identifier>& identifiers, const std::vector<PhaseType>& newPhase, const std::vector<double>& width)>(&MultiInclusions<2>::addLayer))
    .def("changePhase",&MultiInclusions<2>::changePhase)
    .def("getAllIdentifiers",static_cast<vector<Identifier>(MultiInclusions<2>::*)() const>(&MultiInclusions<2>::getAllIdentifiers))
    .def("getAllPhases",static_cast<std::vector<PhaseType>  (MultiInclusions<2>::*)() const>(&MultiInclusions<2>::getAllPhases))
    .def("getAllPhases",static_cast<std::vector<PhaseType> (MultiInclusions<2>::*)(size_t) const>(&MultiInclusions<2>::getAllPhases))
    .def("getAllCenters",&MultiInclusions<2>::getAllCenters)
    .def("getIdentifiers", &MultiInclusions<2>::getIdentifiers)
    .def("setMatrixPhase", &MultiInclusions<2>::setMatrixPhase)
    //.def("translate", &MultiInclusions<2>::translate)
    //.def("setInnerShapes", static_cast<void (MultiInclusions<3>::*)(const vector<smallShape::ConvexPolyhedron<DIM>>&)>(&MultiInclusions<3>::setInnerShapes))
    //.def("setInnerShapes", static_cast<void (MultiInclusions<3>::*)(const vector<smallShape::SphereInc<DIM>>&)>(&MultiInclusions<3>::setInnerShapes))
    //.def("", &MultiInclusions<3>::)
    ;

py::class_<Structure<3>>(merope,"Structure_3D")
	.def(py::init<MultiInclusions<3>>())
	// all same : combination2 initializer
	.def(py::init<Structure<3>, Structure<3>, std::map<PhaseType, PhaseType>>())
    .def(py::init<MultiInclusions<3>, Structure<3>, std::map<PhaseType, PhaseType>>())
    .def(py::init<Structure<3>, MultiInclusions<3>, std::map<PhaseType, PhaseType>>())
    .def(py::init<MultiInclusions<3>, MultiInclusions<3>, std::map<PhaseType, PhaseType>>())
	.def(py::init<Structure<3>, Structure<3>, Structure<3>>())
	.def("getAllPhases", &Structure<3>::getAllPhases)
	;

py::class_<Structure<2>>(merope,"Structure_2D")
    .def(py::init<MultiInclusions<2>>())
    // all same : combination2 initializer
    .def(py::init<Structure<2>, Structure<2>, std::map<PhaseType, PhaseType>>())
    .def(py::init<MultiInclusions<2>, Structure<2>, std::map<PhaseType, PhaseType>>())
    .def(py::init<Structure<2>, MultiInclusions<2>, std::map<PhaseType, PhaseType>>())
    .def(py::init<MultiInclusions<2>, MultiInclusions<2>, std::map<PhaseType, PhaseType>>())
    .def(py::init<Structure<2>, Structure<2>, Structure<2>>())
    .def("getAllPhases", &Structure<2>::getAllPhases)
    ;


py::class_<FieldStructure<3>>(merope, "FieldStructure_3D")
    .def(py::init<const CartesianField<3>&>())
    .def(py::init<const FieldStructure<3>&, const FieldStructure<3>&, const FieldStructure<3>&>())
    .def(py::init<const FieldStructure<3>&, const FieldStructure<3>&, std::function<double(double, double)>>())
    ;

py::class_<FieldStructure<2>>(merope, "FieldStructure_2D")
    .def(py::init<const CartesianField<2>&>())
    .def(py::init<const FieldStructure<2>&, const FieldStructure<2>&, const FieldStructure<2>&>())
    .def(py::init<const FieldStructure<2>&, const FieldStructure<2>&, std::function<double(double, double)>>())
    ;




py::class_<LaguerreTess<3>>(merope,"LaguerreTess_3D")
	.def(py::init<std::array<double,3>, std::vector<Sphere<3>>>())
	.def("setAspRatio", &LaguerreTess<3>::setAspRatio)
    .def("computeTessels", &LaguerreTess<3>::computeTessels)
    .def("getLength",  &LaguerreTess<3>::getL)
	;

py::class_<LaguerreTess<2>>(merope,"LaguerreTess_2D")
    .def(py::init<std::array<double,2>, std::vector<Sphere<2>>>())
    .def("setAspRatio", &LaguerreTess<2>::setAspRatio)
    .def("computeTessels", &LaguerreTess<2>::computeTessels)
    .def("getLength",  &LaguerreTess<2>::getL)
    ;

py::class_<voroInterface::VoroInterface<3>>(merope,"VoroInterface_3D")
	.def(py::init<std::array<double,3>, const std::vector<Sphere<3>>& >())
	.def("drawGnuPlot",&voroInterface::VoroInterface<3>::drawGnuPlot)
	;

py::class_<voroInterface::VoroInterface<2>>(merope,"VoroInterface_2D")
    .def(py::init<std::array<double,2>, const std::vector<Sphere<2>>& >())
    .def("drawGnuPlot",&voroInterface::VoroInterface<2>::drawGnuPlot)
    ;


//merope.def("getClosestNeighbors_2D", &voroInterface::getClosestNeighbors<2>);
merope.def("getClosestNeighbors_3D", &voroInterface::getClosestNeighbors<3>);

//! Gaussian fields

py::module_ gaussianField = merope.def_submodule("gaussianField", "Build (anamorphosed) gaussian field");


py::class_<CartesianField<3>>(merope, "CartesianField_3D")
    //.def(py::init<>())
    .def(py::init<const gaussianField::SimpleGaussianField<3>&, const Point<3>&>())
    .def(py::init<const realScalarField::Field<3>&, const Point<3>&>())
    .def(py::init<const vox::GridField<3>&, const Point<3>&>())
    //.def("set", &CartesianField<3>::setGaussianField)
    //.def("set", &CartesianField<3>::setScalarField)
    //.def("set", &CartesianField<3>::setDiscretizedField)
    .def("getGaussianField", &CartesianField<3>::getGaussianField)
    .def("getScalarField", &CartesianField<3>::getScalarField)
    .def("getDiscretizedField", &CartesianField<3>::getDiscretizedField)
    .def("setLength", &CartesianField<3>::setLength)
    ;

py::class_<CartesianField<2>>(merope, "CartesianField_2D")
    //.def(py::init<>())
    .def(py::init<const gaussianField::SimpleGaussianField<2>&, const Point<2>&>())
    .def(py::init<const realScalarField::Field<2>&, const Point<2>&>())
    .def(py::init<const vox::GridField<2>&, const Point<2>&>())
    //.def("set", &CartesianField<2>::setGaussianField)
    //.def("set", &CartesianField<2>::setScalarField)
    //.def("set", &CartesianField<2>::setDiscretizedField)
    .def("getGaussianField", &CartesianField<2>::getGaussianField)
    .def("getScalarField", &CartesianField<2>::getScalarField)
    .def("getDiscretizedField", &CartesianField<2>::getDiscretizedField)
    .def("setLength", &CartesianField<2>::setLength)
    ;


py::class_<realScalarField::Field<3>>(merope, "ScalarField_3D")
    .def(py::init<std::function<double(Point<3>)>>())
    ;

py::class_<realScalarField::Field<2>>(merope, "ScalarField_2D")
    .def(py::init<std::function<double(Point<2>)>>())
    ;

py::enum_<gaussianField::CovType>(gaussianField,"CovType")
    .value("Spherical",gaussianField::CovType::Spherical)
    .value("Gaussian",gaussianField::CovType::Gaussian)
    .value("Exponential",gaussianField::CovType::Exponential)
    .export_values()
    ;

py::class_<gaussianField::CovSum>(gaussianField, "CovSum")
    .def(py::init<>())
    .def("add", static_cast<void (gaussianField::CovSum::*)(gaussianField::CovType, std::vector<double>)>(&gaussianField::CovSum::add))
    ;

py::class_<gaussianField::StepDis>(gaussianField, "StepDis")
    .def(py::init<>())
    .def(py::init<std::string>())
    ;

py::class_<gaussianField::SimpleGaussianField<3>>(gaussianField, "GaussianField_3D")
    .def(py::init<std::function<double(Point<3>)>, std::function<double(double)>>())
    .def("setNonlinearTransform", &gaussianField::SimpleGaussianField<3>::setNonlinearTransform<gaussianField::StepDis>)
    .def_readwrite("seed", &gaussianField::SimpleGaussianField<3>::seed)
    ;

py::class_<gaussianField::SimpleGaussianField<2>>(gaussianField, "GaussianField_2D")
    .def(py::init<std::function<double(Point<2>)>, std::function<double(double)>>())
    .def("setNonlinearTransform", &gaussianField::SimpleGaussianField<2>::setNonlinearTransform<gaussianField::StepDis>)
    .def_readwrite("seed", &gaussianField::SimpleGaussianField<2>::seed)
    ;


//! Mesh

py::module_ mesh = merope.def_submodule("mesh", "Mesh structure");

py::class_<mesh::generator::MeshGenerator>(mesh, "MeshGenerator")
    .def(py::init<>())
    .def("getMeshOrder",&mesh::generator::MeshGenerator::getMeshOrder)
    .def("setMeshOrder",&mesh::generator::MeshGenerator::setMeshOrder)
    .def("write",static_cast<void(mesh::generator::MeshGenerator::*)(std::string) const>(&mesh::generator::MeshGenerator::write))
    .def("setMultiInclusions",&mesh::generator::MeshGenerator::setMultiInclusions)
    .def("getMeshSize",&mesh::generator::MeshGenerator::getMeshSize)
    .def("setMeshSize",&mesh::generator::MeshGenerator::setMeshSize)
    .def("setAdimMergeDistance0",&mesh::generator::MeshGenerator::setAdimMergeDistance0)
    .def("setAdimMergeDistance1",&mesh::generator::MeshGenerator::setAdimMergeDistance1)
    .def("setBinaryOutput", &mesh::generator::MeshGenerator::setBinaryOutput)
    ;

//! MicroInclusions

py::module_ microInclusion = merope.def_submodule("microInclusion", "Functions for micro inclusions");

py::class_<merope::smallShape::ConvexPolyhedronInc<3>>(microInclusion, "ConvexPolyhedronInc_3D")
    ;

py::class_<merope::smallShape::ConvexPolyhedronInc<2>>(microInclusion, "ConvexPolyhedronInc_2D")
    ;

py::class_<merope::smallShape::SpheroPolyhedronInc<3>>(microInclusion, "SpheroPolyhedronInc_3D")
    ;

py::class_<merope::smallShape::SpheroPolyhedronInc<2>>(microInclusion, "SpheroPolyhedronInc_2D")
    ;

py::class_<smallShape::PolyhedronFactory<3>>(microInclusion,"PolyhedronFactory_3D")
    .def(py::init<>())
    .def("fromVertices", & smallShape::PolyhedronFactory<3>::fromVertices)
    ;

py::class_<smallShape::PolyhedronFactory<2>>(microInclusion,"PolyhedronFactory_2D")
    .def(py::init<>())
    .def("fromVertices", & smallShape::PolyhedronFactory<2>::fromVertices)
    ;

py::class_<smallShape::SpheroPolyhedronFactory<3>>(microInclusion,"SpheroPolyhedronFactory_3D")
    .def(py::init<>())
    .def("fromVertices", & smallShape::SpheroPolyhedronFactory<3>::fromVertices)
    ;

py::class_<smallShape::SpheroPolyhedronFactory<2>>(microInclusion,"SpheroPolyhedronFactory_2D")
    .def(py::init<>())
    .def("fromVertices", & smallShape::SpheroPolyhedronFactory<2>::fromVertices)
    ;

//! Simplified interfaces
py::class_<smallShape::LayerInstructions>(merope,"LayerInstructions")
     .def(py::init<Identifier, PhaseType, double>())
     .def_readwrite("identifier", &smallShape::LayerInstructions::identifier)
     .def_readwrite("phase", &smallShape::LayerInstructions::phase)
     .def_readwrite("width", &smallShape::LayerInstructions::width)
     ;

py::class_<Pre_InterfaceMultiInclusions<3>>(merope,"Pre_InterfaceMultiInclusions_3D")
    // here only for SimpleStructure_3D. Do not allow for much methods
    .def(py::init<>())
    .def("setSpheres",static_cast<void (Pre_InterfaceMultiInclusions<3>::*)(std::vector<Sphere<3>>)>(&Pre_InterfaceMultiInclusions<3>::setSpheres))
    .def("getSpheres",&Pre_InterfaceMultiInclusions<3>::getSpheres)
    .def("fromFile",static_cast<void (Pre_InterfaceMultiInclusions<3>::*)(std::string)>(&Pre_InterfaceMultiInclusions<3>::fromFile))
    .def("setAspRatio",&Pre_InterfaceMultiInclusions<3>::setAspRatio)
    .def("setTypeCrystal",&Pre_InterfaceMultiInclusions<3>::setTypeCrystal)
    ;

py::class_<Pre_InterfaceMultiInclusions<2>>(merope,"Pre_InterfaceMultiInclusions_2D")
    // here only for SimpleStructure_3D. Do not allow for much methods
    .def(py::init<>())
    .def("setSpheres",static_cast<void (Pre_InterfaceMultiInclusions<2>::*)(std::vector<Sphere<2>>)>(&Pre_InterfaceMultiInclusions<2>::setSpheres))
    .def("getSpheres",&Pre_InterfaceMultiInclusions<2>::getSpheres)
    .def("fromFile",static_cast<void (Pre_InterfaceMultiInclusions<2>::*)(std::string)>(&Pre_InterfaceMultiInclusions<2>::fromFile))
    .def("setAspRatio",&Pre_InterfaceMultiInclusions<2>::setAspRatio)
    .def("setTypeCrystal",&Pre_InterfaceMultiInclusions<2>::setTypeCrystal)
    ;

py::class_<InterfaceMultiInclusions<3>,Pre_InterfaceMultiInclusions<3>>(merope,"SimpleMultiInclusions_3D")
    .def(py::init<>())
    .def("setLayerList", &InterfaceMultiInclusions<3>::setLayerList)
    .def("getLayerList", &InterfaceMultiInclusions<3>::getLayerList)
    .def("build", &InterfaceMultiInclusions<3>::build)
    .def("setLength", &InterfaceMultiInclusions<3>::setLength)
    ;

py::class_<InterfaceMultiInclusions<2>,Pre_InterfaceMultiInclusions<2>>(merope,"SimpleMultiInclusions_2D")
    .def(py::init<>())
    .def("setLayerList", &InterfaceMultiInclusions<2>::setLayerList)
    .def("getLayerList", &InterfaceMultiInclusions<2>::getLayerList)
    .def("build", &InterfaceMultiInclusions<2>::build)
    .def("setLength", &InterfaceMultiInclusions<2>::setLength)
    ;


py::class_<InterfaceStructure<3>>(merope,"SimpleStructure_3D")
   .def(py::init<>())
   .def("setErosionPhase", &InterfaceStructure<3>::setErosionPhase)
   .def("setErosionInclusionsPhase", &InterfaceStructure<3>::setErosionInclusionsPhase)
   .def("setErosionWidth", &InterfaceStructure<3>::setErosionWidth)
   .def("setInnerPhase", &InterfaceStructure<3>::setInnerPhase)
   .def("setColorization", &InterfaceStructure<3>::setColorization)
   .def("setLength", &InterfaceStructure<3>::setLength)
   .def("build", &InterfaceStructure<3>::build)
   .def("frac2erosionWidth", &InterfaceStructure<3>::frac2erosionWidth)
   .def("erosionWidth2frac", &InterfaceStructure<3>::erosionWidth2frac)
   .def_readwrite("mainInclusions",&InterfaceStructure<3>::mainInclusions)
   .def_readwrite("secdInclusions",&InterfaceStructure<3>::secdInclusions)
   ;

py::class_<InterfaceStructure<2>>(merope,"SimpleStructure_2D")
   .def(py::init<>())
   .def("setErosionPhase", &InterfaceStructure<2>::setErosionPhase)
   .def("setErosionInclusionsPhase", &InterfaceStructure<2>::setErosionInclusionsPhase)
   .def("setErosionWidth", &InterfaceStructure<2>::setErosionWidth)
   .def("setInnerPhase", &InterfaceStructure<2>::setInnerPhase)
   .def("setColorization", &InterfaceStructure<2>::setColorization)
   .def("setLength", &InterfaceStructure<2>::setLength)
   .def("build", &InterfaceStructure<2>::build)
   .def("frac2erosionWidth", &InterfaceStructure<2>::frac2erosionWidth)
   .def("erosionWidth2frac", &InterfaceStructure<2>::erosionWidth2frac)
   .def_readwrite("mainInclusions",&InterfaceStructure<2>::mainInclusions)
   .def_readwrite("secdInclusions",&InterfaceStructure<2>::secdInclusions)
   ;

//!

merope.def("printRandomOrient_3D", geneOrientation::printRandomOrient_3D);

//! Enums

py::enum_<TypeCrystal>(merope,"TypeCrystal")
    .value("Voronoi",TypeCrystal::Voronoi)
    .value("Laguerre",TypeCrystal::Laguerre)
    .value("JohnsonMehl",TypeCrystal::JohnsonMehl)
    .value("Spheres",TypeCrystal::Spheres)
    .export_values()
    ;

py::enum_<ColorMaterialID>(merope, "ColorMaterialID")
    .value("Poly",ColorMaterialID::Poly)
    .value("Erode",ColorMaterialID::Erode)
    .value("Erode2Mat",ColorMaterialID::Erode2Mat)
    .value("Erode3Mat",ColorMaterialID::Erode3Mat)
    .export_values()
    ;

py::enum_<homogenization::Rule>(merope,"HomogenizationRule")
	.value("Reuss", homogenization::Rule::Reuss)
	.value("Voigt", homogenization::Rule::Voigt)
	.value("Smallest", homogenization::Rule::Smallest)
	.value("Largest", homogenization::Rule::Largest)
	.export_values()
	;

py::enum_<vox::VoxelRule>(merope,"VoxelRule")
	.value("Average", 	vox::VoxelRule::Average)
	.value("Center", 	vox::VoxelRule::Center)
	.export_values()
	;

//! fonctions
merope.def("setNbOfThreads",merope::setNbOfThreads);
merope.def("symmetrize", static_cast<void(*)(std::string, std::string, std::array<size_t, 3>)>(&vox::symmetrize<3>));
merope.def("symmetrize", static_cast<void(*)(std::string, std::string, std::array<size_t, 2>)>(&vox::symmetrize<2>));
}

PYBIND11_MODULE(Merope, merope_){
	createModule_merope(merope_);
}

PYBIND11_MODULE(merope, merope_){
	createModule_merope(merope_);
}

