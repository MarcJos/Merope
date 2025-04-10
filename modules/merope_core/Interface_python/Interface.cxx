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
#include <functional>

#include "../../AlgoPacking/include/AlgoNames.hxx"
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
#include "../include/Voxellation/DynamicVoxellizer.hxx"
#include "../include/VTKinout/VTK_printer_for_dynamic_voxellizer.hxx"
#include "../include/Voxellation/GridAnalyzer.hxx"
#include "../include/Voxellation/VoxRecurStructure.hxx"
#include "../include/Voxellation/Deprecated/Voxellation.hxx"
#include "../include/geneOrientations.hxx"
#include "../include/VTKinout/VTKRead.hxx"
#include "../include/VTKinout/TIFFRead.hxx"
#include "../include/VTKinout/VTKStream.hxx"
#include "../include/Grid/GridTypes.hxx"
#include "../include/MicroToMesh/MeshGenerator.hxx"
#include "../include/MicroInclusion/MicroInclusion.hxx"
#include "../include/AlgoLaguerre/MakeCentroidal.hxx"
#include "../include/AlgoLaguerre/Optimize_LaguerreTess.hxx"

#include "Interf_function_pointer.hxx"

namespace py = pybind11;

using namespace merope;
using namespace sac_de_billes;


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
        1, std::vector<std::size_t> {myVector.size()}, std::vector<std::size_t> {sizeof(T)}));
}

template<class VOX>
py::array_t<double> computeField_Numpy(VOX& vox) {
    vector<double> result{};
    {// stop the GIL of python just for this computation
        Stop_gil stop_gil{};
        result = vox.getField();
    }
    return create_matrix(result);
}

template<unsigned short DIM>
class NumpyConverter {
public:
    py::array_t<double> compute_RealField(const vox::voxellizer::GridRepresentation<DIM>& input) {
        return create_matrix(input.template get<double>());
    }
    py::array_t<sac_de_billes::PhaseType> compute_PhaseField(const vox::voxellizer::GridRepresentation<DIM>& input) {
        return create_matrix(input.template get<sac_de_billes::PhaseType>());
    }
};

template<class T>
py::list process_list(const std::vector<T>& vecteur) {
    py::list output_list;
    for (const auto& item : vecteur) {
        output_list.append(item);
    }
    return output_list;
}

template<unsigned short DIM>
py::list get_as_list(const vox::voxellizer::GridRepresentation<DIM>& grid_rep) {
    auto copy_of_grid_rep = grid_rep;
    copy_of_grid_rep.convert_to_stl_format();
    py::list result{};
    copy_of_grid_rep.try_on_all([&](const auto& compo) {
        result = process_list(compo);
        });
    return result;
}

template<class VOX>
py::array_t<merope::PhaseType> computePurePhaseGrid(VOX& vox) {
    vector<merope::PhaseType> result{};
    {// stop the GIL of python just for this computation
        Stop_gil stop_gil{};
        std::vector<double> coeffs{};
        result = vox.computePurePhaseGrid(coeffs);
    }
    return create_matrix(result);
}


template<unsigned short DIM>
inline void create_functions_depending_on_dimension(py::module_& merope, std::string DIM_S) {

    py::module_ geometry = merope.def_submodule("geometry", "Basic geometry primitives and functions");
    py::module_ vox = merope.def_submodule("vox", "Functions for voxellation");
    py::module_ microInclusion = merope.def_submodule("microInclusion", "Functions for micro inclusions");
    py::module_ gaussianField = merope.def_submodule("gaussianField", "Build (anamorphosed) gaussian field");


    auto my_string_Segment = std::string("Segment") + DIM_S;
    py::class_<sac_de_billes::Segment<DIM>>(geometry, my_string_Segment.c_str())
        .def(py::init<array<Point<DIM>, 2>>())
        ;

    auto my_string_HalfSpace = std::string("HalfSpace") + DIM_S;
    py::class_<merope::HalfSpace<DIM>>(geometry, my_string_HalfSpace.c_str())
        .def(py::init<Point<DIM>, double>())
        .def(py::init<Point<DIM>, Point<DIM>>())
        .def("c", static_cast<double(merope::HalfSpace<DIM>::*)()const>(&merope::HalfSpace<DIM>::c))
        .def("vec", &merope::HalfSpace<DIM>::vec)
        ;


    auto my_string_GridField = "GridField" + DIM_S;
    py::class_<merope::vox::GridField<DIM>>(vox, my_string_GridField.c_str())
        ;

    auto my_string_get_linear_index = "get_linear_index" + DIM_S;
    vox.def(my_string_get_linear_index.c_str(), &merope::vox::auxi::get_linear_index<DIM, array<std::size_t, DIM>>);


    vox.def("symmetrize", static_cast<void(*)(std::string, std::string, std::array<std::size_t, DIM>)>(&merope::vox::symmetrize<DIM>));

    auto my_string_CartesianField = "CartesianField" + DIM_S;
    py::class_<CartesianField<DIM>>(vox, my_string_CartesianField.c_str())
        //.def(py::init<>())
        .def(py::init<const gaussianField::SimpleGaussianField<DIM>&, const Point<DIM>&>())
        .def(py::init<const gaussianField::NumericalCovariance<DIM>&, const Point<DIM>&>())
        .def(py::init<const realScalarField::Field<DIM>&, const Point<DIM>&>())
        .def(py::init<const merope::vox::GridField<DIM>&, const Point<DIM>&>())
        //.def("set", &CartesianField<DIM>::setGaussianField)
        //.def("set", &CartesianField<DIM>::setScalarField)
        //.def("set", &CartesianField<DIM>::setDiscretizedField)
        .def("getGaussianField", &CartesianField<DIM>::getGaussianField)
        .def("getScalarField", &CartesianField<DIM>::getScalarField)
        .def("getDiscretizedField", &CartesianField<DIM>::getDiscretizedField)
        .def("setLength", &CartesianField<DIM>::setLength)
        ;

    auto my_string_vtk_printer = "vtk_printer" + DIM_S;
    py::class_<vtk_adapter::vtk_printer<DIM>>(vox, my_string_vtk_printer.c_str())
        .def(py::init<>())
        .def("printVTK", &vtk_adapter::vtk_printer<DIM>::printVTK,
            py::arg("gridRepresentation"), py::arg("fileVTK"), py::arg("nameValue") = "MaterialId")
        .def("printVTK_segmented", &vtk_adapter::vtk_printer<DIM>::printVTK_segmented,
            py::arg("gridRepresentation"), py::arg("fileVTK"), py::arg("fileCoeff"), py::arg("nameValue") = "MaterialId")
        .def("printVTK_removeUnusedPhase", &vtk_adapter::vtk_printer<DIM>::printVTK_removeUnusedPhase,
            py::arg("gridRepresentation"), py::arg("fileVTK"), py::arg("fileCoeff"), py::arg("nameValue") = "MaterialId")
        ;

    auto my_string_GridRepresentation = "GridRepresentation" + DIM_S;
    py::class_<vox::voxellizer::GridRepresentation<DIM>>(vox, my_string_GridRepresentation.c_str())
        .def(py::init<const Structure<DIM>&, vox::GridParameters<DIM>, vox::VoxelRule, std::map<tuple<PhaseType, PhaseType>, PhaseType>>())
        .def(py::init<const Structure<DIM>&, vox::GridParameters<DIM>, vox::VoxelRule>())
        .def(py::init<const FieldStructure<DIM>&, vox::GridParameters<DIM>, vox::VoxelRule>())
        .def(py::init<const vox::voxellizer::GridRepresentation<DIM>&>()) // deep copy
        .def("apply_coefficients", &vox::voxellizer::GridRepresentation<DIM>::apply_coefficients)
        .def("apply_homogRule", static_cast<void(vox::voxellizer::GridRepresentation<DIM>::*)(homogenization::Rule)> (&vox::voxellizer::GridRepresentation<DIM>::apply_homogRule))
        .def("apply_homogRule", static_cast<void(vox::voxellizer::GridRepresentation<DIM>::*)(homogenization::Rule, std::vector<double>)> (&vox::voxellizer::GridRepresentation<DIM>::apply_homogRule))
        .def("apply_texture", [](vox::voxellizer::GridRepresentation<DIM>& gridRep, Interf_TexturePointer i1) {
        gridRep.apply_texture(i1.get<DIM, 1>());
            })
        .def("convert_to_stl_format", &vox::voxellizer::GridRepresentation<DIM>::convert_to_stl_format)
        .def("convert_to_Iso_format", &vox::voxellizer::GridRepresentation<DIM>::convert_to_Iso_format)
        .def("convert_to_AnIso_format", &vox::voxellizer::GridRepresentation<DIM>::convert_to_AnIso_format)
        .def("removeUnusedPhase", &vox::voxellizer::GridRepresentation<DIM>::removeUnusedPhase)
        .def("get_PureRealField", &vox::voxellizer::GridRepresentation<DIM>::template get<vox::composite::Pure<double>>)
        .def("get_PurePhaseField", &vox::voxellizer::GridRepresentation<DIM>::template get<vox::composite::Pure<PhaseType>>)
        .def("get_as_list", &get_as_list<DIM>)
        .def("__str__", &vox::voxellizer::GridRepresentation<DIM>::to_string_type)
        ;

    auto my_string_GridAnalyzer = "GridAnalyzer" + DIM_S;
    py::class_<vox::grid_analyzer::GridAnalyzer<DIM>>(vox, my_string_GridAnalyzer.c_str())
        .def(py::init<>())
        .def("print_percentages", &vox::grid_analyzer::GridAnalyzer<DIM>::print_percentages)
        .def("compute_percentages", &vox::grid_analyzer::GridAnalyzer<DIM>::compute_percentages)
        ;

    auto my_string_GridParameters = "GridParameters" + DIM_S;
    py::class_<vox::GridParameters<DIM>>(vox, my_string_GridParameters.c_str())
        .def(py::init<std::array<size_t, DIM>, std::array<double, DIM>>())
        ;

    auto my_string_create_grid_parameters_N_L = "create_grid_parameters_N_L" + DIM_S;
    vox.def(my_string_create_grid_parameters_N_L.c_str(), static_cast<vox::GridParameters<DIM>(*)(std::array<size_t, DIM>, std::array<double, DIM>)> (vox::create_grid_parameters_N_L<DIM>));
    auto my_string_create_sub_grid_parameters_N_L = "create_grid_parameters_N_L" + DIM_S;
    vox.def(my_string_create_sub_grid_parameters_N_L.c_str(), static_cast<vox::GridParameters<DIM>(*)(std::array<size_t, DIM>, std::array<size_t, DIM>, std::array<size_t, DIM>, std::array<double, DIM>)> (vox::create_grid_parameters_N_L<DIM>));


    auto my_string_NumpyConverter = "NumpyConverter" + DIM_S;
    py::class_<NumpyConverter<DIM>>(vox, my_string_NumpyConverter.c_str())
        .def(py::init<>())
        .def("compute_RealField", &NumpyConverter<DIM>::compute_RealField)
        .def("compute_PhaseField", &NumpyConverter<DIM>::compute_PhaseField)
        ;

    auto my_string_Voxellation = std::string("Voxellation") + DIM_S;
    py::class_<vox::Voxellation<DIM>>(merope, my_string_Voxellation.c_str())
        .def(py::init<const MultiInclusions<DIM>&>())
        .def(py::init<const Structure<DIM>&>())
        .def(py::init<const FieldStructure<DIM>&>())
        .def("proceed", static_cast<void(vox::Voxellation<DIM>::*)(std::array<std::size_t, DIM>, std::array<std::size_t, DIM>, std::array<std::size_t, DIM>)>(&vox::Voxellation<DIM>::proceed), py::call_guard<Stop_gil>())
        .def("proceed", static_cast<void(vox::Voxellation<DIM>::*)(std::array<std::size_t, DIM>)>(&vox::Voxellation<DIM>::proceed), py::call_guard<Stop_gil>())
        .def("computePhaseGrid", &vox::Voxellation<DIM>::computePhaseGrid, py::call_guard<Stop_gil>())
        .def("computeCompositeGrid", &vox::Voxellation<DIM>::computeCompositeGrid, py::call_guard<Stop_gil>())
        .def("printFile", &vox::Voxellation<DIM>::printFile, py::call_guard<Stop_gil>())
        .def("getField_Numpy", &computeField_Numpy<vox::Voxellation<DIM>>)
        .def("getField", &vox::Voxellation<DIM>::getField)
        .def("getPhaseGrid", &computePurePhaseGrid<vox::Voxellation<DIM>>)
        .def("setPureCoeffs", &vox::Voxellation<DIM>::setPureCoeffs)
        .def("setHomogRule", &vox::Voxellation<DIM>::setHomogRule)
        .def("setVoxelRule", &vox::Voxellation<DIM>::setVoxelRule)
        .def("printFieldFile", &vox::Voxellation<DIM>::printFieldFile)
        .def("__str__", &vox::Voxellation<DIM>::print)
        ;

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
        .def("enlarge", static_cast<void (MultiInclusions<DIM>::*)(const std::vector<Identifier>&identifiers, double width)>(&MultiInclusions<DIM>::enlarge))
        .def("enlarge", static_cast<void (MultiInclusions<DIM>::*)(const std::vector<Identifier>&identifiers, const std::vector<double>&width)>(&MultiInclusions<DIM>::enlarge))
        .def("addLayer", static_cast<void (MultiInclusions<DIM>::*)(const std::vector<Identifier>&identifiers, PhaseType newPhase, double width)>(&MultiInclusions<DIM>::addLayer))
        .def("addLayer", static_cast<void (MultiInclusions<DIM>::*)(const std::vector<Identifier>&identifiers, const std::vector<PhaseType>&newPhase, const std::vector<double>&width)>(&MultiInclusions<DIM>::addLayer))
        .def("changePhase", &MultiInclusions<DIM>::changePhase)
        .def("getAllIdentifiers", static_cast<vector<Identifier>(MultiInclusions<DIM>::*)() const>(&MultiInclusions<DIM>::getAllIdentifiers))
        .def("getAllPhases", static_cast<std::vector<PhaseType>(MultiInclusions<DIM>::*)() const>(&MultiInclusions<DIM>::getAllPhases))
        .def("getAllPhases", static_cast<std::vector<PhaseType>(MultiInclusions<DIM>::*)(std::size_t) const>(&MultiInclusions<DIM>::getAllPhases))
        .def("getAllCenters", &MultiInclusions<DIM>::getAllCenters)
        .def("getIdentifiers", &MultiInclusions<DIM>::getIdentifiers)
        .def("setMatrixPhase", &MultiInclusions<DIM>::setMatrixPhase)
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
        .def(py::init(
            [](const FieldStructure<DIM>& a1, const FieldStructure<DIM>& a2, Interf_FuncPointer i1) {
                auto func = i1.get<2, 1>();
                std::function<double(double, double)> func2 = [func](double x1, double x2) {
                    Point<2> pt = { x1, x2 };
                    return func(pt);
                    };
                return FieldStructure<DIM>(a1, a2, func2);
            }
        ))
        .def(py::init(
            [](const FieldStructure<DIM>& a1, const FieldStructure<DIM>& a2, std::function<double(double, double)> a3) {
                error_use_omp();
                return FieldStructure<DIM>(a1, a2, a3);
            }))
        ;

    auto my_string_LaguerreTess = "LaguerreTess" + DIM_S;
    py::class_<LaguerreTess<DIM>>(merope, my_string_LaguerreTess.c_str())
        .def(py::init<std::array<double, DIM>, std::vector<Sphere<DIM>>>())
        .def("setAspRatio", &LaguerreTess<DIM>::setAspRatio)
        .def("computeTessels", &LaguerreTess<DIM>::computeTessels)
        .def("getLength", &LaguerreTess<DIM>::getL)
        .def("toPolyInclusions", &LaguerreTess<DIM>::toPolyInclusions)
        ;

    auto my_string_VoroInterface = "VoroInterface" + DIM_S;
    py::class_<voroInterface::VoroInterface<DIM>>(merope, my_string_VoroInterface.c_str())
        .def(py::init<std::array<double, DIM>, const std::vector<Sphere<DIM>>& >())
        .def(py::init<std::array<double, DIM>, const std::vector<Sphere<DIM>>&, std::array<bool, DIM> >())
        .def("findTessel", &voroInterface::VoroInterface<DIM>::findTessel)
        .def("computeSolids", &voroInterface::VoroInterface<DIM>::computeSolids)
        .def("drawGnuPlot", &voroInterface::VoroInterface<DIM>::drawGnuPlot)
        .def("drawCellsPov", &voroInterface::VoroInterface<DIM>::drawCellsPov)
        .def("printCustom", &voroInterface::VoroInterface<DIM>::printCustom)
        .def("addWallCylinder", &voroInterface::VoroInterface<DIM>::addWallCylinder)
        .def("getCellCenters", &voroInterface::VoroInterface<DIM>::getCellCenters)
        ;

    auto my_string_ScalarField = "ScalarField" + DIM_S;
    py::class_<realScalarField::Field<DIM>>(merope, my_string_ScalarField.c_str())
        .def(py::init(
            [](Interf_FuncPointer i1) {
                return realScalarField::Field<DIM>(i1.get<DIM, 1>());
            }
        ))
        .def(py::init(
            [](std::function<double(Point<DIM>)> f) {
                error_use_omp();
                return realScalarField::Field<DIM>(f);
            }))
        ;


    auto my_string_GaussianField = "GaussianField" + DIM_S;
    py::class_<gaussianField::SimpleGaussianField<DIM>>(gaussianField, my_string_GaussianField.c_str())
        .def(py::init(
            [](Interf_FuncPointer i1, Interf_FuncPointer i2) {
                return gaussianField::SimpleGaussianField<DIM>(i1.get<DIM, 1>(), i2.get<1, 1>());
            }
        ))
        .def(py::init(
            [](std::function<double(Point<DIM>)> covariance_,
                std::function<double(double)> nonlinearFunction_) {
                    error_use_omp();
                    return gaussianField::SimpleGaussianField<DIM>(covariance_, nonlinearFunction_);
            }
        ))
        .def("setNonlinearTransform", &gaussianField::SimpleGaussianField<DIM>::template setNonlinearTransform<gaussianField::StepDis>)
        .def_readwrite("seed", &gaussianField::SimpleGaussianField<DIM>::seed)
        ;

    auto my_string_NumericalCovariance = "NumericalCovariance" + DIM_S;
    py::class_<gaussianField::NumericalCovariance<DIM>>(gaussianField, my_string_NumericalCovariance.c_str())
        .def(py::init(
            [](Interf_FuncPointer i1) {
                return gaussianField::NumericalCovariance<DIM>(i1.get<DIM, 1>());
            }
        ))
        .def(py::init(
            [](std::function<double(Point<DIM>)> f) {
                error_use_omp();
                return gaussianField::NumericalCovariance<DIM>(f);
            }
        ))
        ;



}


inline void createModule_merope(py::module_& merope) {

    //! Gaussian fields
    py::module_ geometry = merope.def_submodule("geometry", "Basic geometry primitives and functions");
    py::module_ vox = merope.def_submodule("vox", "Functions for voxellation");
    py::module_ gaussianField = merope.def_submodule("gaussianField", "Build (anamorphosed) gaussian field");
    py::module_ mesh = merope.def_submodule("mesh", "Mesh structure");

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
        .def("phaseFraction", &Grid_VER::phaseFraction)
        .def("__str__", static_cast<void(Grid_VER::*)() const>(&Grid_VER::print))
        ;

    //merope.def("getClosestNeighbors_2D", &voroInterface::getClosestNeighbors<2>);
    merope.def("getClosestNeighbors_3D", &voroInterface::getClosestNeighbors<3>);


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

    py::enum_<mesh::gmsh_writer::MeshMethod>(mesh, "MeshMethod")
        .value("native_gmsh", mesh::gmsh_writer::MeshMethod::native_gmsh)
        .value("OpenCascade", mesh::gmsh_writer::MeshMethod::OpenCascade)
        .export_values()
        ;

    py::class_<mesh::generator::MeshGenerator>(mesh, "MeshGenerator")
        .def(py::init<>())
        .def("getMeshOrder", &mesh::generator::MeshGenerator::getMeshOrder)
        .def("setMeshOrder", &mesh::generator::MeshGenerator::setMeshOrder)
        .def("write", static_cast<void(mesh::generator::MeshGenerator::*)(std::string, mesh::gmsh_writer::MeshMethod) const>(&mesh::generator::MeshGenerator::write),
            py::arg("nameFile") = "myFile.geo", py::arg("meshMethod") = mesh::gmsh_writer::MeshMethod::native_gmsh)
        .def("setMultiInclusions", &mesh::generator::MeshGenerator::setMultiInclusions)
        .def("getMeshSize", &mesh::generator::MeshGenerator::getMeshSize)
        .def("setMeshSize", &mesh::generator::MeshGenerator::setMeshSize)
        .def("setAdimMergeDistance", &mesh::generator::MeshGenerator::setAdimMergeDistance)
        .def("setBinaryOutput", &mesh::generator::MeshGenerator::setBinaryOutput)
        .def("do_not_mesh", &mesh::generator::MeshGenerator::do_not_mesh)
        .def("set_nameOutput", &mesh::generator::MeshGenerator::set_nameOutput)
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

    py::enum_<merope::vox::VoxelRule>(vox, "VoxelRule")
        .value("Average", merope::vox::VoxelRule::Average)
        .value("Center", merope::vox::VoxelRule::Center)
        .value("Laminate", merope::vox::VoxelRule::Laminate)
        .value("PolyGeom", merope::vox::VoxelRule::PolyGeom)
        .export_values()
        ;

    //! fonctions
    merope.def("setNbOfThreads", omp_set_num_threads);
    merope.def("getNbOfThreads", sac_de_billes::auxi_function::get_num_threads);

    merope.def("makeCentroidal_3D", merope::optimizeLaguerreTess::makeCentroidal<3, 3>);

    py::class_<Interf_FuncPointer>(merope, "Interf_FuncPointer")
        .def(py::init<std::size_t, std::array<std::size_t, 2>>())
        ;

    py::class_<Interf_TexturePointer>(merope, "Interf_TexturePointer")
        .def(py::init<std::size_t, std::array<std::size_t, 2>>())
        ;

    py::class_<sac_de_billes::Cylinder<3>>(geometry, "Cylinder")
        .def(py::init<const Segment<3>&, double>())
        .def("volume", &sac_de_billes::Cylinder<3>::volume)
        .def_readwrite("radius", &sac_de_billes::Cylinder<3>::radius)
        .def_readwrite("axis", &sac_de_billes::Cylinder<3>::axis)
        ;

    py::class_<CylinderInclusions<3>>(merope, "CylinderInclusions_3D")
        .def(py::init<>())
        .def("setLength", &CylinderInclusions<3>::setLength)
        .def("addInclusion", &CylinderInclusions<3>::addInclusion)
        .def("setInclusions", &CylinderInclusions<3>::setInclusions)
        ;

}


inline void create_last_function(py::module_& merope) {

    py::module_ microInclusion = merope.def_submodule("microInclusion", "Functions for micro inclusions");

    py::class_<merope::smallShape::CylinderInc<3>>(microInclusion, "CylinderInclusions")
        .def(py::init<const Cylinder<3> &>())
        ;

    // Adding a new method afterwards
    py::class_<MultiInclusions<3>> MultiInclusions_3D_binding = merope.attr("MultiInclusions_3D");
    MultiInclusions_3D_binding.def("setInclusions", static_cast<void (MultiInclusions<3>::*)(const CylinderInclusions<3>&)>(&MultiInclusions<3>::setInclusions));
}

PYBIND11_MODULE(Merope, merope_) {
    omp_set_num_threads(4);
    createModule_merope(merope_);
    create_functions_depending_on_dimension<3>(merope_, std::string("_3D"));
    create_functions_depending_on_dimension<2>(merope_, std::string("_2D"));
    create_last_function(merope_);
}

PYBIND11_MODULE(merope, merope_) {
    omp_set_num_threads(4);
    createModule_merope(merope_);
    create_functions_depending_on_dimension<3>(merope_, std::string("_3D"));
    create_functions_depending_on_dimension<2>(merope_, std::string("_2D"));
    create_last_function(merope_);
}



