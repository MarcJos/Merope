//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../Grid/CartesianGrid.hxx"
#include "../Grid/VtkHistogramm.hxx"



namespace merope {

namespace vox {
namespace convertGrid {


//! \return a grid converting each voxel according to the transformation rule
//! \param  grid0 : input grid
//! \param rule : a function for converting each voxel
template<bool use_open_mp, unsigned short DIM, class C1, class C2, class FUNCTION>
CartesianGrid<DIM, C1> localConvert(const CartesianGrid<DIM, C2>& grid0, const FUNCTION& rule);

//! @return a grid converting each voxel according to the texturer
//! @tparam TEXTURER : a function of the voxel AND the spatial position x
//! @param grid0 : input grid
//! @param texturer : a function (Point<DIM> x, C2 voxelphase)
template<unsigned short DIM, class C1, class C2, class TEXTURER>
CartesianGrid<DIM, C1> applyFunctionDependingOnX(const CartesianGrid<DIM, C2>& grid0, const TEXTURER& texturer);

//! transform a grid according to the transformation rule
//! \param : grid0, input grid
//! \param rule : a function applied on each voxel
template<unsigned short DIM, class C2, class FUNCTION>
void apply_inplace(CartesianGrid<DIM, C2>& grid0, FUNCTION rule);

//! natural conversion of all voxels
template<unsigned short DIM, class COMPOSITE_OUT, class COMPOSITE_IN>
CartesianGrid<DIM, COMPOSITE_OUT> convert_to(const CartesianGrid<DIM, COMPOSITE_IN>& grid0);
//! \return a grid of volume fractions from a phase grid
//! \param grid0 : grid to be converted
template<unsigned short DIM, class COMPOSITE>
CartesianGrid<DIM, composite::to_stl_format<COMPOSITE>> convert_to_stl_format(const CartesianGrid<DIM, COMPOSITE>& grid0);

//! \return linearize a Cartesian Grid
//! \param grid0 : grid to be converted
template<unsigned short DIM, class VOXEL_TYPE>
vector<VOXEL_TYPE> linearize(const CartesianGrid<DIM, VOXEL_TYPE>& grid0);

//! \return a grid of (unsigned short) phases grom a field, by building an order histogramm of it
//! \param grid0 : grid to be converted
//! \param fieldValues : track the values
template<unsigned short DIM>
CartesianGrid<DIM, PhaseType> fromFieldToPhase(const CartesianGrid<DIM, double>& gridField, vector<double>& fieldValues);

//! takes a grid of phases with associated coefficients
//! changes the values of the phases and fieldValues such that the grid0 contains phases from 0->N and that fieldValues is increasing
template<unsigned short DIM, class COEFF_TYPE, class = enable_if_t<is_arithmetic_v<COEFF_TYPE>>>
void renormalizeWithCoefficients(CartesianGrid<DIM, PhaseType>& grid0, vector<COEFF_TYPE>& fieldValues);

//! @brief : transform a given grid of PhaseType of values inside [[0, N]] with possibly zero voxel of phase
//! n \in [[0, N]] into a grid of PhaseType of values inside [[0, N']], with, for each n in [[0, N']]
//! there exists at least one voxel of phase n
//! phase order is preserved
//! @param grid : inplace modified grid
template<unsigned short DIM>
void removeUnusedPhase(CartesianGrid<DIM, PhaseType>& grid, vector<PhaseType>& phaseValues);

//! slightly less that 16000, constant for avoiding overflow in CartesianGrid<DIM, unsigned short>;
static constexpr unsigned short NB_PHASE_USHORT = 15995;

namespace auxi {
//! \return a grid converting each voxel according to the transformation rule
//! \param : grid0, input grid
//! \param rule : a function for converting each voxel
template<bool use_open_mp, unsigned short DIM, class C1, class C2, class FUNCTION>
CartesianGrid<DIM, C1> localConvertCartesian(const CartesianGrid<DIM, C2>& grid0, FUNCTION rule);
}  // namespace auxi

}  // namespace convertGrid
}  // namespace vox
}  // namespace merope

#include "../Grid/ConvertGrix.ixx"


