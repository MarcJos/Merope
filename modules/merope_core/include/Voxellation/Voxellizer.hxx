//! Copyright : see license.txt
//!
//! \briefClass for voxellization of inclusions
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Geometry/GeomTools.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../Grid/ConvertGrix.hxx"
#include "../Grid/GridManipulations.hxx"
#include "../Grid/GridTypes.hxx"
#include "../Grid/GridTypesBase.hxx"
#include "../Grid/VtkHistogramm.hxx"
#include "../Voxellation/VoxRecurStructure.hxx"
#include "../VTKinout/VTKStream.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<class STRUCTURE>
using Phase_Type_From_Structure_Type = std::conditional_t<std::is_same_v<STRUCTURE, Structure<2>> or std::is_same_v<STRUCTURE, Structure<3>>, PhaseType, double>;

namespace voxellizer {
//! retrieve the natural output of a Structure + VoxelRule
template<unsigned short DIM, VoxelRule VOXEL_RULE, class STRUCTURE>
vox::CartesianGrid<DIM, composite::OutputFormat<VOXEL_RULE, DIM, Phase_Type_From_Structure_Type<STRUCTURE>>> transformStructIntoGrid(const STRUCTURE& structure, GridParameters<DIM> preSubGrid) {
    return transformIntoGrid<DIM, composite::OutputFormat<VOXEL_RULE, DIM, Phase_Type_From_Structure_Type<STRUCTURE>>>(structure, preSubGrid);
}

//! @return : the result of the application of a homogenization rule
//! @param phaseFracVol : a cartesian grid of composite isotropic voxel
template<homogenization::Rule HOMOG_RULE, unsigned short DIM>
GridField<DIM> applyHomogRule_T(const CartesianGrid<DIM, composite::Iso<double>>& phaseFracVol);
//! @return : the result of the application of a homogenization rule
//! @param phaseFracVol : a cartesian grid of composite isotropic voxel
//! @param coeffs : coefficients
template<homogenization::Rule HOMOG_RULE, unsigned short DIM>
GridField<DIM> applyHomogRule_T(const CartesianGrid<DIM, composite::Iso<PhaseType>>& phaseFracVol,
    const vector<double>& pureCoeffs);

//! @return : the result of an application of a texturer value <- f(x, phase)
//! @param grid0 : initial grid of composite phases
//! @param texturer :  a function f : \R^DIM x PHASE_ENSEMBLE -> \R
template<unsigned short DIM, class COMPOSITE, class TEXTURER, typename = enable_if_t<composite::is_composite<COMPOSITE>>>
CartesianGrid<DIM, composite::Change_Type_Composite<DIM, COMPOSITE, double>> apply_texture(const CartesianGrid<DIM, COMPOSITE>& grid0, const TEXTURER& texturer);

//! @return : the grid with applied coefficients from a composite grid of phase
//! @param grid : grid assumed to be of phase
//! @param pureCoeffs : coefficient to transform the phase into a real
template<unsigned short DIM, class COMPOSITE>
CartesianGrid<DIM, composite::Change_Type_Composite<DIM, COMPOSITE, double>> apply_coefficients(const CartesianGrid<DIM, COMPOSITE>& grid, const vector<double>& pureCoeffs);

}  // namespace  voxellizer

} /* namespace vox */
}  // namespace merope

#include "../Voxellation/Voxellizer.ixx"