//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRIDMANIPULATIONS_HXX_
#define GRIDMANIPULATIONS_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"


#include "../Grid/GridTypes.hxx"
#include "../Grid/VoxelRule.hxx"
#include "../Grid/InclusionAndVoxel.hxx"
#include "../Grid/ConvertGrix.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

namespace vox {
namespace gridAuxi {

//! \return an array of {volumeFraction, coefficients}
//! \param phaseFrac : a list of phase and volume fractions
//! \param coeff :  a list of coefficients corresponding the the phases
//! \fixme : Il vaudrait mieux traiter d'une part phaseFrac, et considérer coeff ailleurs.
array<vector<double>, 2> getTabCoeff(const VoxelPhaseFrac& phaseFrac,
        const vector<double>& coeff);
//! \return the combination of two grids by means of a function
//! \param grid1, grid2 : voxellations
//! \param function : function depending on elements of grid1 and grid2
//! \see vox::gridAuxi::combineVoxelFunc
template<unsigned short DIM, class VOXEL_TYPE, class FUNCTION>
CartesianGrid<DIM, VOXEL_TYPE> combineGridFunction(const CartesianGrid<DIM, VOXEL_TYPE>& grid1,
        const CartesianGrid<DIM, VOXEL_TYPE>& grid2,
        const FUNCTION  func);
//! \return the combination of two grids by means of a mask
//! \param grid1, grid2 : voxellations
//! \param mask : if 1 triggers replacement of grid1 by grid2
//! \see vox::gridAuxi::combineVoxelMask
template<unsigned short DIM, class VOXEL_TYPE>
CartesianGrid<DIM, VOXEL_TYPE> combineGridMask(const CartesianGrid<DIM, VOXEL_TYPE>& grid1,
        const CartesianGrid<DIM, VOXEL_TYPE>& grid2,
        const CartesianGrid<DIM, VOXEL_TYPE>& mask);
//! \return the result on a single voxel of a mask operation
//! \param voxelData1, voxelData2 : voxel data to be masked
//! \param voxelMask : each phase = 1 triggers replacement of voxelData1 by voxelData2
template<class VOXEL_TYPE>
VOXEL_TYPE combineVoxelMask(const VOXEL_TYPE& voxelData1, const VOXEL_TYPE& voxelData2, const VOXEL_TYPE& voxelMask);
//! \return the result on a single voxel of a function operation
//! \param voxelData1, voxelData2 : voxel data to be combined
//! \param function : the resulting phase(s) are given by function(phase1, phase2)
template<class VOXEL_TYPE, class FUNCTION>
VOXEL_TYPE combineVoxelFunc(const VOXEL_TYPE& voxelData1, const VOXEL_TYPE& voxelData2, const FUNCTION& function);
//! \return the coefficient of the convex combination to be applied on phases densities
//! \param voxelMask : data of the voxel in the mask
//! \seed combineVoxelMask
array<double, 2> translateMask(const VoxelPhaseFrac& voxelMask);

} // namespace gridAuxi
} // namespace vox
} // namespace merope

#include "../Grid/GridManipulations.ixx"

#endif /* GRIDMANIPULATIONS_HXX_ */
