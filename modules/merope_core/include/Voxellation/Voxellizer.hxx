//! Copyright : see license.txt
//!
//! \briefClass for voxellization of inclusions
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/GeomTools.hxx"

#include "../Grid/CartesianGrid.hxx"
#include "../Grid/ConvertGrix.hxx"
#include "../Grid/GridManipulations.hxx"
#include "../Grid/GridTypes.hxx"
#include "../SingleVoxel/SingleVoxel_Headers.hxx"
#include "../Grid/VtkHistogramm.hxx"
#include "../Voxellation/VoxRecurStructure.hxx"
#include "../VTKinout/VTKStream.hxx"


namespace merope {
namespace vox {

namespace voxellizer {
template<unsigned short DIM, class VOXEL_POLICY, class STRUCTURE>
vox::CartesianGrid<DIM, composite::OutputFormat<VOXEL_POLICY::voxelRule, DIM, Phase_Type_From_Structure_Type<STRUCTURE>>>
transformStructIntoGrid(const STRUCTURE& structure, GridParameters<DIM> preSubGrid, VOXEL_POLICY voxelPolicy) {
    return transformIntoGrid<DIM>(structure, preSubGrid, voxelPolicy);
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