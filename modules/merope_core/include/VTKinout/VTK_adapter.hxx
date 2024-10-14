//! Copyright : see license.txt
//!
//! \brief

#pragma once


#include "../MeropeNamespace.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "VTKStream.hxx"

namespace merope {
namespace vtk_adapter {

//! @return name of vtk type corresponding to PHASE_OUT
template<class PHASE_OUT>
std::string name_type();

template<class PHASE_OUT, unsigned short DIM, class PHASE_IN>
void printVTK_removeUnusedPhase(const vox::CartesianGrid<DIM, PHASE_IN>& cartesianGrid,
    string fileVTK,
    string fileCoeffs,
    string nameValue = "MaterialId");

template<class PHASE_OUT, unsigned short DIM, class PHASE_IN>
void printVTK(const vox::CartesianGrid<DIM, PHASE_IN>& cartesianGrid,
    string fileVTK,
    string nameValue = "MaterialId");

template<class PHASE_OUT, unsigned short DIM, class PHASE_IN>
void printVTK_segmented(const vox::CartesianGrid<DIM, PHASE_IN>& cartesianGrid,
    string fileVTK, string fileCoeff,
    string nameValue = "MaterialId");

template<class PHASE_OUT, unsigned short DIM, class PHASE_IN>
void printVTK(VTKstream& vtkstream,
    array<size_t, DIM> n, array<double, DIM> dx,
    const vector<PHASE_IN>& values, string nameValue);

template<unsigned short DIM, class PHASE_COEFFS>
void printCoeffs(const vector<PHASE_COEFFS>& coefficients, string fileCoeff);

}  // namespace  vtk_adapter

}  // namespace  merope

#include "VTK_adapter.ixx"

