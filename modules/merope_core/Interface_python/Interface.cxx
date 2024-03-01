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
#include "../include/AlgoLaguerre/MakeCentroidal.hxx"
#include "../include/AlgoLaguerre/Optimize_LaguerreTess.hxx"

namespace py = pybind11;

using namespace merope;
using namespace sac_de_billes;

// Voxellation


// stop the GIL of Python in order to allow Python functions to be used in C++
// beware, this is dangerous if Python creates/deletes objects
class Stop_gil {
public:
    Stop_gil() :_save{ nullptr } {
        _save = PyEval_SaveThread();
    }
    ~Stop_gil() {
        PyEval_RestoreThread(_save);
    }

private:
    PyThreadState* _save;
};

template<class T>
py::array_t<T> create_matrix(vector<T> myVector) {
    return py::array_t<T>(py::buffer_info(myVector.data(), sizeof(T), py::format_descriptor<T>::format(),
        1, std::vector<size_t> {myVector.size()}, std::vector<size_t> {sizeof(T)}));
}

template<class VOX>
py::array_t<double> getField_Numpy(VOX& vox) {
    vector<double> result{};
    {// stop the GIL of python just for this computation
        Stop_gil stop_gil{};
        result = vox.getField();
    }
    return create_matrix(result);
}

template<class VOX>
py::array_t<vox::VTK_PHASE> getPhaseGrid(VOX& vox) {
    vector<vox::VTK_PHASE> result{};
    {// stop the GIL of python just for this computation
        Stop_gil stop_gil{};
        result = vox.getPhaseGrid();
    }
    return create_matrix(result);
}


template<unsigned short DIM>
inline void create_functions_depending_on_dimension(py::module_& merope, std::string DIM_S) {

    auto my_string_Voxellation = std::string("Voxellation") + DIM_S;
    py::class_<vox::Voxellation<DIM>>(merope, my_string_Voxellation.c_str())
        .def(py::init<const MultiInclusions<DIM>&>())
        .def(py::init<const Structure<DIM>&>())
        .def(py::init<const FieldStructure<DIM>&>())
        .def("proceed", static_cast<void(vox::Voxellation<DIM>::*)(std::array<size_t, DIM>, std::array<size_t, DIM>, std::array<size_t, DIM>)>(&vox::Voxellation<DIM>::proceed), py::call_guard<Stop_gil>())
        .def("proceed", static_cast<void(vox::Voxellation<DIM>::*)(std::array<size_t, DIM>)>(&vox::Voxellation<DIM>::proceed), py::call_guard<Stop_gil>())
        .def("computePhaseGrid", &vox::Voxellation<DIM>::computePhaseGrid, py::call_guard<Stop_gil>())
        .def("printFile", &vox::Voxellation<DIM>::printFile, py::call_guard<Stop_gil>())
        .def("getField_Numpy", &getField_Numpy<vox::Voxellation<DIM>>)
        .def("getField", &vox::Voxellation<DIM>::getField)
        .def("getPhaseGrid", &getPhaseGrid<vox::Voxellation<DIM>>)
        .def("setPureCoeffs", &vox::Voxellation<DIM>::setPureCoeffs)
        .def("setHomogRule", &vox::Voxellation<DIM>::setHomogRule)
        .def("setVoxelRule", &vox::Voxellation<DIM>::setVoxelRule)
        .def("printFieldFile", &vox::Voxellation<DIM>::printFieldFile)
        ;

    auto my_string_GridField = "GridField" + DIM_S;
    py::class_<vox::GridField<DIM>>(merope, my_string_GridField.c_str())
        ;

    auto my_string_get_linear_index = "get_linear_index" + DIM_S;
    merope.def(my_string_get_linear_index.c_str(), &vox::auxi::get_linear_index<DIM, array<size_t, DIM>>);

    //! Microstructure
    auto my_string_SphereInclusions = "SphereInclusions" + DIM_S;
    py::class_<SphereInclusions<DIM>>(merope, my_string_SphereInclusions.c_str())
        .def(py::init<>())
        .def("setSpheres", static_cast<void(SphereInclusions<DIM>::*)(const std::vector<Sphere<DIM>>)>(&SphereInclusions<DIM>::setSpheres))
        .def("fromFile", static_cast<void(SphereInclusions<DIM>::*)(std::string)>(&SphereInclusions<DIM>::fromFile))
        .def("fromFile", static_cast<void(SphereInclusions<DIM>::*)(std::string, std::string)>(&SphereInclusions<DIM>::fromFile))
        .def("setLength", &SphereInclusions<DIM>::setLength)
        .def("getLength", &SphereInclusions<DIM>::getLength)
        .def("printDump", &SphereInclusions<DIM>::printDump)
        .def("fromHisto", &SphereInclusions<DIM>::fromHisto)
        .def("printFracVol", &SphereInclusions<DIM>::printFracVol)
        .def("printDump", &SphereInclusions<DIM>::printDump)
        .def("printPos", &SphereInclusions<DIM>::printPos)
        .def("printCSV", &SphereInclusions<DIM>::printCSV_space)
        .def("printVER", &SphereInclusions<DIM>::printVER)
        .def("getSpheres", &SphereInclusions<DIM>::getSpheres)
        .def("translate", &SphereInclusions<DIM>::translate)
        ;

    auto my_string_PolyInclusions = "PolyInclusions" + DIM_S;
    py::class_<PolyInclusions<DIM>>(merope, my_string_PolyInclusions.c_str())
        .def(py::init<>())
        .def("setLength", &PolyInclusions<DIM>::setLength)
        .def("addInclusion", &PolyInclusions<DIM>::addInclusion)
        .def("setInclusions", &PolyInclusions<DIM>::setInclusions)
        ;

    //! MicroInclusions
    py::module_ microInclusion = merope.def_submodule("microInclusion", "Functions for micro inclusions");

    auto my_string_ConvexPolyhedronInc = "ConvexPolyhedronInc" + DIM_S;
    py::class_<merope::smallShape::ConvexPolyhedronInc<DIM>>(microInclusion, my_string_ConvexPolyhedronInc.c_str())
        ;

    auto my_string_SpheroPolyhedronInc = "SpheroPolyhedronInc" + DIM_S;
    py::class_<merope::smallShape::SpheroPolyhedronInc<DIM>>(microInclusion, my_string_SpheroPolyhedronInc.c_str())
        ;

    auto my_string_PolyhedronFactory = "PolyhedronFactory" + DIM_S;
    py::class_<smallShape::PolyhedronFactory<DIM>>(microInclusion, my_string_PolyhedronFactory.c_str())
        .def(py::init<>())
        .def("fromVertices", &smallShape::PolyhedronFactory<DIM>::fromVertices)
        ;

    auto my_string_SpheroPolyhedronFactory = "SpheroPolyhedronFactory" + DIM_S;
    py::class_<smallShape::SpheroPolyhedronFactory<DIM>>(microInclusion, my_string_SpheroPolyhedronFactory.c_str())
        .def(py::init<>())
        .def("fromVertices", &smallShape::SpheroPolyhedronFactory<DIM>::fromVertices)
        ;

    auto my_string_Pre_InterfaceMultiInclusions = "Pre_InterfaceMultiInclusions" + DIM_S;
    py::class_<Pre_InterfaceMultiInclusions<DIM>>(merope, my_string_Pre_InterfaceMultiInclusions.c_str())
        // here only for SimpleStructure_DIMD. Do not allow for much methods
        .def(py::init<>())
        .def("setSpheres", static_cast<void (Pre_InterfaceMultiInclusions<DIM>::*)(std::vector<Sphere<DIM>>)>(&Pre_InterfaceMultiInclusions<DIM>::setSpheres))
        .def("getSpheres", &Pre_InterfaceMultiInclusions<DIM>::getSpheres)
        .def("fromFile", static_cast<void (Pre_InterfaceMultiInclusions<DIM>::*)(std::string)>(&Pre_InterfaceMultiInclusions<DIM>::fromFile))
        .def("setAspRatio", &Pre_InterfaceMultiInclusions<DIM>::setAspRatio)
        .def("setTypeCrystal", &Pre_InterfaceMultiInclusions<DIM>::setTypeCrystal)
        ;

    auto my_string_SimpleMultiInclusions = "SimpleMultiInclusions" + DIM_S;
    py::class_<InterfaceMultiInclusions<DIM>, Pre_InterfaceMultiInclusions<DIM>>(merope, my_string_SimpleMultiInclusions.c_str())
        .def(py::init<>())
        .def("setLayerList", &InterfaceMultiInclusions<DIM>::setLayerList)
        .def("getLayerList", &InterfaceMultiInclusions<DIM>::getLayerList)
        .def("build", &InterfaceMultiInclusions<DIM>::build)
        .def("setLength", &InterfaceMultiInclusions<DIM>::setLength)
        ;

    auto my_string_SimpleStructure = "SimpleStructure" + DIM_S;
    py::class_<InterfaceStructure<DIM>>(merope, my_string_SimpleStructure.c_str())
        .def(py::init<>())
        .def("setErosionPhase", &InterfaceStructure<DIM>::setErosionPhase)
        .def("setErosionInclusionsPhase", &InterfaceStructure<DIM>::setErosionInclusionsPhase)
        .def("setErosionWidth", &InterfaceStructure<DIM>::setErosionWidth)
        .def("setInnerPhase", &InterfaceStructure<DIM>::setInnerPhase)
        .def("setColorization", &InterfaceStructure<DIM>::setColorization)
        .def("setLength", &InterfaceStructure<DIM>::setLength)
        .def("build", &InterfaceStructure<DIM>::build)
        .def("frac2erosionWidth", &InterfaceStructure<DIM>::frac2erosionWidth)
        .def("erosionWidth2frac", &InterfaceStructure<DIM>::erosionWidth2frac)
        .def_readwrite("mainInclusions", &InterfaceStructure<DIM>::mainInclusions)
        .def_readwrite("secdInclusions", &InterfaceStructure<DIM>::secdInclusions)
        ;

    auto my_string_algo_fit_volumes = "algo_fit_volumes" + DIM_S;
    py::class_<merope::optimizeLaguerreTess::algo_fit_volumes<DIM>>(merope, my_string_algo_fit_volumes.c_str())
        .def(py::init<const Point<DIM>&, const vector<Sphere<DIM>>&, const vector<double>&>())
        .def("proceed", &merope::optimizeLaguerreTess::algo_fit_volumes<DIM>::proceed)
        .def("getCenterTessels", &merope::optimizeLaguerreTess::algo_fit_volumes<DIM>::getCenterTessels)
        .def("maxDeltaVolumes", &merope::optimizeLaguerreTess::algo_fit_volumes<DIM>::maxDeltaVolumes)
        .def("getCurrentVolumes", &merope::optimizeLaguerreTess::algo_fit_volumes<DIM>::getCurrentVolumes)
        ;

    auto my_string_SpheroPolyInclusions = "SpheroPolyInclusions" + DIM_S;
    py::class_<SpheroPolyhedronInclusions<DIM>>(merope, my_string_SpheroPolyInclusions.c_str())
        .def(py::init<>())
        .def("setLength", &SpheroPolyhedronInclusions<DIM>::setLength)
        .def("addInclusion", &SpheroPolyhedronInclusions<DIM>::addInclusion)
        .def("setInclusions", &SpheroPolyhedronInclusions<DIM>::setInclusions)
        ;

    auto my_string_MultiInclusions = "MultiInclusions" + DIM_S;
    py::class_<MultiInclusions<DIM>>(merope, my_string_MultiInclusions.c_str())
        .def(py::init<>())
        .def("setInclusions", static_cast<void (MultiInclusions<DIM>::*)(LaguerreTess<DIM>)>(&MultiInclusions<DIM>::setInclusions))
        .def("setInclusions", static_cast<void (MultiInclusions<DIM>::*)(const PolyInclusions<DIM>&)>(&MultiInclusions<DIM>::setInclusions))
        .def("setInclusions", static_cast<void (MultiInclusions<DIM>::*)(const SphereInclusions<DIM>&)>(&MultiInclusions<DIM>::setInclusions))
        .def("setInclusions", static_cast<void (MultiInclusions<DIM>::*)(const EllipseInclusions<DIM>&)>(&MultiInclusions<DIM>::setInclusions))
        .def("setInclusions", static_cast<void (MultiInclusions<DIM>::*)(const std::vector<smallShape::SpheroPolyhedronInc<DIM>>&, Point<DIM>)>(&MultiInclusions<DIM>::setInclusions))
        .def("setInclusions", static_cast<void (MultiInclusions<DIM>::*)(const SpheroPolyhedronInclusions<DIM>&)>(&MultiInclusions<DIM>::setInclusions))
        .def("addLayer", static_cast<void (MultiInclusions<DIM>::*)(const std::vector<Identifier>&identifiers, PhaseType newPhase, double width)>(&MultiInclusions<DIM>::addLayer))
        .def("addLayer", static_cast<void (MultiInclusions<DIM>::*)(const std::vector<Identifier>&identifiers, const std::vector<PhaseType>&newPhase, const std::vector<double>&width)>(&MultiInclusions<DIM>::addLayer))
        .def("changePhase", &MultiInclusions<DIM>::changePhase)
        .def("getAllIdentifiers", static_cast<vector<Identifier>(MultiInclusions<DIM>::*)() const>(&MultiInclusions<DIM>::getAllIdentifiers))
        .def("getAllPhases", static_cast<std::vector<PhaseType>(MultiInclusions<DIM>::*)() const>(&MultiInclusions<DIM>::getAllPhases))
        .def("getAllPhases", static_cast<std::vector<PhaseType>(MultiInclusions<DIM>::*)(size_t) const>(&MultiInclusions<DIM>::getAllPhases))
        .def("getAllCenters", &MultiInclusions<DIM>::getAllCenters)
        .def("getIdentifiers", &MultiInclusions<DIM>::getIdentifiers)
        .def("setMatrixPhase", &MultiInclusions<DIM>::setMatrixPhase)
        //.def("translate", &MultiInclusions<DIM>::translate)
        //.def("setInnerShapes", static_cast<void (MultiInclusions<DIM>::*)(const vector<smallShape::ConvexPolyhedron<DIM>>&)>(&MultiInclusions<DIM>::setInnerShapes))
        //.def("setInnerShapes", static_cast<void (MultiInclusions<DIM>::*)(const vector<smallShape::SphereInc<DIM>>&)>(&MultiInclusions<DIM>::setInnerShapes))
        //.def("", &MultiInclusions<DIM>::)
        ;

    auto my_string_Structure = "Structure" + DIM_S;
    py::class_<Structure<DIM>>(merope, my_string_Structure.c_str())
        .def(py::init<MultiInclusions<DIM>>())
        // all same : combination2 initializer
        .def(py::init<Structure<DIM>, Structure<DIM>, std::map<PhaseType, PhaseType>>())
        .def(py::init<MultiInclusions<DIM>, Structure<DIM>, std::map<PhaseType, PhaseType>>())
        .def(py::init<Structure<DIM>, MultiInclusions<DIM>, std::map<PhaseType, PhaseType>>())
        .def(py::init<MultiInclusions<DIM>, MultiInclusions<DIM>, std::map<PhaseType, PhaseType>>())
        .def(py::init<Structure<DIM>, Structure<DIM>, Structure<DIM>>())
        .def("getAllPhases", &Structure<DIM>::getAllPhases)
        ;

    auto my_string_FieldStructure = "FieldStructure" + DIM_S;
    py::class_<FieldStructure<DIM>>(merope, my_string_FieldStructure.c_str())
        .def(py::init<const CartesianField<DIM>&>())
        .def(py::init<const FieldStructure<DIM>&, const FieldStructure<DIM>&, const FieldStructure<DIM>&>())
        .def(py::init<const FieldStructure<DIM>&, const FieldStructure<DIM>&, std::function<double(double, double)>>())
        ;

    auto my_string_LaguerreTess = "LaguerreTess" + DIM_S;
    py::class_<LaguerreTess<DIM>>(merope, my_string_LaguerreTess.c_str())
        .def(py::init<std::array<double, DIM>, std::vector<Sphere<DIM>>>())
        .def("setAspRatio", &LaguerreTess<DIM>::setAspRatio)
        .def("computeTessels", &LaguerreTess<DIM>::computeTessels)
        .def("getLength", &LaguerreTess<DIM>::getL)
        ;

    auto my_string_VoroInterface = "VoroInterface" + DIM_S;
    py::class_<voroInterface::VoroInterface<DIM>>(merope, my_string_VoroInterface.c_str())
        .def(py::init<std::array<double, DIM>, const std::vector<Sphere<DIM>>& >())
        .def("drawGnuPlot", &voroInterface::VoroInterface<DIM>::drawGnuPlot)
        ;

    py::module_ gaussianField = merope.def_submodule("gaussianField", "Build (anamorphosed) gaussian field");

    auto my_string_CartesianField = "CartesianField" + DIM_S;
    py::class_<CartesianField<DIM>>(merope, my_string_CartesianField.c_str())
        //.def(py::init<>())
        .def(py::init<const gaussianField::SimpleGaussianField<DIM>&, const Point<DIM>&>())
        .def(py::init<const gaussianField::NumericalCovariance<DIM>&, const Point<DIM>&>())
        .def(py::init<const realScalarField::Field<DIM>&, const Point<DIM>&>())
        .def(py::init<const vox::GridField<DIM>&, const Point<DIM>&>())
        //.def("set", &CartesianField<DIM>::setGaussianField)
        //.def("set", &CartesianField<DIM>::setScalarField)
        //.def("set", &CartesianField<DIM>::setDiscretizedField)
        .def("getGaussianField", &CartesianField<DIM>::getGaussianField)
        .def("getScalarField", &CartesianField<DIM>::getScalarField)
        .def("getDiscretizedField", &CartesianField<DIM>::getDiscretizedField)
        .def("setLength", &CartesianField<DIM>::setLength)
        ;

    auto my_string_ScalarField = "ScalarField" + DIM_S;
    py::class_<realScalarField::Field<DIM>>(merope, my_string_ScalarField.c_str())
        .def(py::init<std::function<double(Point<DIM>)>>())
        ;

    auto my_string_GaussianField = "GaussianField" + DIM_S;
    py::class_<gaussianField::SimpleGaussianField<DIM>>(gaussianField, my_string_GaussianField.c_str())
        .def(py::init<std::function<double(Point<DIM>)>, std::function<double(double)>>())
        .def("setNonlinearTransform", &gaussianField::SimpleGaussianField<DIM>::template setNonlinearTransform<gaussianField::StepDis>)
        .def_readwrite("seed", &gaussianField::SimpleGaussianField<DIM>::seed)
        ;

    auto my_string_NumericalCovariance = "NumericalCovariance" + DIM_S;
    py::class_<gaussianField::NumericalCovariance<DIM>>(gaussianField, my_string_NumericalCovariance.c_str())
        .def(py::init<std::function<double(Point<DIM>)>>())
        ;

    merope.def("symmetrize", static_cast<void(*)(std::string, std::string, std::array<size_t, DIM>)>(&vox::symmetrize<DIM>));


}


void createModule_merope(py::module_& merope) {
    //VTKinout
    py::class_<VTKstream>(merope, "VTKstream")
        .def(py::init<const char*>())
        .def("close", &VTKstream::close)
        ;

    py::class_<VTKRead>(merope, "VTKRead")
        .def(py::init<const char*>())
        ;

    py::class_<TIFFRead>(merope, "TIFFRead")
        .def(py::init<const char*>())
        ;

    py::class_<Grid_VER>(merope, "Grid")
        .def(py::init<>())
        .def(py::init<VTKRead&>())
        .def("DivNormFactor", &Grid_VER::DivNormFactor)
        .def("toVTKCELL", &Grid_VER::toVTKCELL)
        .def("VTKheaderCELL", &Grid_VER::VTKheaderCELL)
        .def("phaseFraction", &Grid_VER::phaseFraction)
        ;

    //merope.def("getClosestNeighbors_2D", &voroInterface::getClosestNeighbors<2>);
    merope.def("getClosestNeighbors_3D", &voroInterface::getClosestNeighbors<3>);

    //! Gaussian fields
    py::module_ gaussianField = merope.def_submodule("gaussianField", "Build (anamorphosed) gaussian field");

    py::enum_<gaussianField::CovType>(gaussianField, "CovType")
        .value("Spherical", gaussianField::CovType::Spherical)
        .value("Gaussian", gaussianField::CovType::Gaussian)
        .value("Exponential", gaussianField::CovType::Exponential)
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

    //! Mesh

    py::module_ mesh = merope.def_submodule("mesh", "Mesh structure");

    py::class_<mesh::generator::MeshGenerator>(mesh, "MeshGenerator")
        .def(py::init<>())
        .def("getMeshOrder", &mesh::generator::MeshGenerator::getMeshOrder)
        .def("setMeshOrder", &mesh::generator::MeshGenerator::setMeshOrder)
        .def("write", static_cast<void(mesh::generator::MeshGenerator::*)(std::string) const>(&mesh::generator::MeshGenerator::write))
        .def("setMultiInclusions", &mesh::generator::MeshGenerator::setMultiInclusions)
        .def("getMeshSize", &mesh::generator::MeshGenerator::getMeshSize)
        .def("setMeshSize", &mesh::generator::MeshGenerator::setMeshSize)
        .def("setAdimMergeDistance0", &mesh::generator::MeshGenerator::setAdimMergeDistance0)
        .def("setAdimMergeDistance1", &mesh::generator::MeshGenerator::setAdimMergeDistance1)
        .def("setBinaryOutput", &mesh::generator::MeshGenerator::setBinaryOutput)
        ;


    //! Simplified interfaces
    py::class_<smallShape::LayerInstructions>(merope, "LayerInstructions")
        .def(py::init<Identifier, PhaseType, double>())
        .def_readwrite("identifier", &smallShape::LayerInstructions::identifier)
        .def_readwrite("phase", &smallShape::LayerInstructions::phase)
        .def_readwrite("width", &smallShape::LayerInstructions::width)
        ;
    //!

    merope.def("printRandomOrient_3D", geneOrientation::printRandomOrient_3D);

    //! Enums

    py::enum_<TypeCrystal>(merope, "TypeCrystal")
        .value("Voronoi", TypeCrystal::Voronoi)
        .value("Laguerre", TypeCrystal::Laguerre)
        .value("JohnsonMehl", TypeCrystal::JohnsonMehl)
        .value("Spheres", TypeCrystal::Spheres)
        .export_values()
        ;

    py::enum_<ColorMaterialID>(merope, "ColorMaterialID")
        .value("Poly", ColorMaterialID::Poly)
        .value("Erode", ColorMaterialID::Erode)
        .value("Erode2Mat", ColorMaterialID::Erode2Mat)
        .value("Erode3Mat", ColorMaterialID::Erode3Mat)
        .export_values()
        ;

    py::enum_<homogenization::Rule>(merope, "HomogenizationRule")
        .value("Reuss", homogenization::Rule::Reuss)
        .value("Voigt", homogenization::Rule::Voigt)
        .value("Smallest", homogenization::Rule::Smallest)
        .value("Largest", homogenization::Rule::Largest)
        .export_values()
        ;

    py::enum_<vox::VoxelRule>(merope, "VoxelRule")
        .value("Average", vox::VoxelRule::Average)
        .value("Center", vox::VoxelRule::Center)
        .export_values()
        ;

    //! fonctions
    merope.def("setNbOfThreads", merope::setNbOfThreads);
    merope.def("makeCentroidal_3D", merope::optimizeLaguerreTess::makeCentroidal<3, 3>);
}

PYBIND11_MODULE(Merope, merope_) {
    createModule_merope(merope_);
    create_functions_depending_on_dimension<3>(merope_, std::string("_3D"));
    create_functions_depending_on_dimension<2>(merope_, std::string("_2D"));
}

PYBIND11_MODULE(merope, merope_) {
    createModule_merope(merope_);
    create_functions_depending_on_dimension<3>(merope_, std::string("_3D"));
    create_functions_depending_on_dimension<2>(merope_, std::string("_2D"));
}

