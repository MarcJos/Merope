//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../Grid/ConvertGrix.hxx"

namespace merope {

namespace vtk_adapter {

template<class PHASE_OUT>
std::string name_type() {
    if constexpr (std::is_same_v<PHASE_OUT, unsigned short>) {
        return "unsigned_short";
    }
    if constexpr (std::is_same_v<PHASE_OUT, float>) {
        return "float";
    }
    if constexpr (std::is_same_v<PHASE_OUT, double>) {
        return "double";
    }
    throw std::invalid_argument("PHASE_OUT");
}

template<class PHASE_OUT, unsigned short DIM, class PHASE_IN>
void printVTK_removeUnusedPhase(const vox::CartesianGrid<DIM, PHASE_IN>& cartesianGrid,
    string fileVTK,
    string fileCoeffs,
    string nameValue) {
    auto copy_cartesianGrid = cartesianGrid;
    vector<PhaseType> coefficients = {};
    vox::convertGrid::removeUnusedPhase(copy_cartesianGrid, coefficients);
    printVTK<PHASE_OUT>(copy_cartesianGrid, fileVTK, nameValue);
    printCoeffs<DIM>(coefficients, fileCoeffs);
}

template<class PHASE_OUT, unsigned short DIM, class PHASE_IN>
void printVTK(const vox::CartesianGrid<DIM, PHASE_IN>& cartesianGrid, string fileVTK,
    string nameValue) {
    auto n = cartesianGrid.getNbNodeSubGrid();
    auto dx = cartesianGrid.getDx();
    VTKstream vtkstream(fileVTK.c_str());
    printVTK<PHASE_OUT, DIM, PHASE_IN>(vtkstream, n, dx, cartesianGrid, nameValue);
}

template<class PHASE_OUT, unsigned short DIM, class PHASE_IN>
void printVTK(VTKstream& vtkstream,
    array<size_t, DIM> n, array<double, DIM> dx,
    const vector<PHASE_IN>& values, string nameValue) {
    // header
    vtkstream.STRUCTURED_POINTS<DIM>(n, dx);
    //
    vtkstream.setCELL(values.size());
    //
    if constexpr (std::is_arithmetic_v<PHASE_IN>) {
        vtkstream << "SCALARS " << nameValue << " " << name_type<PHASE_OUT>() << endl;
        vtkstream << "LOOKUP_TABLE default" << endl;
        vtkstream.writeVector<PHASE_OUT>(values, n);
    } else if constexpr (is_same_v<PHASE_IN, array<double, DIM>>) {
        vtkstream << "FIELD FieldData 1" << endl;
        vtkstream << nameValue << " " << DIM << " " << values.size() << " " << name_type<PHASE_OUT>() << endl;
        vtkstream.writeVector<PHASE_OUT>(values, n);
    }
}

template<class PHASE_OUT, unsigned short DIM, class PHASE_IN>
void printVTK_segmented(const vox::CartesianGrid<DIM, PHASE_IN>& cartesianGrid,
    string fileVTK, string fileCoeff, string nameValue) {
    vector<double> coefficients{};
    auto pure_grid = vox::convertGrid::fromFieldToPhase(cartesianGrid, coefficients);
    vox::convertGrid::renormalizeWithCoefficients(pure_grid, coefficients);
    printVTK<PHASE_OUT>(pure_grid, fileVTK, nameValue);
    printCoeffs<DIM>(coefficients, fileCoeff);
}

template<unsigned short DIM, class PHASE_COEFFS>
void printCoeffs(const vector<PHASE_COEFFS>& coefficients, string fileCoeff) {
    ofstream ost(fileCoeff);
    for (auto c : coefficients) {
        ost << c << endl;
    }
}

}  // namespace  vtk_adapter

}  // namespace  merope

