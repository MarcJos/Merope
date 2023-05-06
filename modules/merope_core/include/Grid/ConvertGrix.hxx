//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_CONVERTGRIX_HXX_
#define GRID_CONVERTGRIX_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Grid/CartesianGrid.hxx"
#include "../Grid/VtkHistogramm.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

namespace vox {
namespace convertGrid {


//! \return a grid converting each voxel according to the transformation rule
//! \param : grid0, input grid
//! \param rule : a function for converting each voxel
template<unsigned short DIM, class C1, class C2, class FUNCTION>
CartesianGrid<DIM, C1> localConvert(const CartesianGrid<DIM, C2>& grid0, const FUNCTION& rule);

//! transform a grid according to the transformation rule
//! \param : grid0, input grid
//! \param rule : a function applied on each voxel
template<unsigned short DIM, class C2, class FUNCTION>
void apply_inplace(CartesianGrid<DIM, C2>& grid0, FUNCTION rule);


//! \return a grid of volume fractions from a phase grid
//! \param grid0 : grid to be converted
template<unsigned short DIM>
CartesianGrid<DIM, VoxelPhaseFrac> fromPhaseToFracVol(const CartesianGrid<DIM, VTK_PHASE>& grid0);

//! \return a grid of volume fractions from a phase grid
//! \param grid0 : grid to be converted
template<unsigned short DIM>
vector<vector<tuple<vox::VTK_PHASE, double>>> fromCartesianToVector(CartesianGrid<DIM, VoxelPhaseFrac> grid0);

//! \return a grid of (unsigned short) phases grom a field, by building an order histogramm of it
//! \param grid0 : grid to be converted
//! \param fieldValues : track the values
template<unsigned short DIM>
CartesianGrid<DIM, VTK_PHASE> fromFieldToPhase(const CartesianGrid<DIM, double>& gridField, vector<double>& fieldValues);

//! takes a grid of phases with associated coefficients
//! changes the values of the phases and fieldValues such that the grid0 contains phases from 0->N and that fieldValues is increasing
template<unsigned short DIM>
void renormalizeWithCoefficients(CartesianGrid<DIM, VTK_PHASE>& grid0, vector<double>& fieldValues);

//! slightly less that 16000, constant for avoiding overflow in CartesianGrid<DIM, unsigned short>;
static constexpr unsigned short NB_PHASE_USHORT = 15995;

namespace auxi {
//! \return a grid converting each voxel according to the transformation rule
//! \param : grid0, input grid
//! \param rule : a function for converting each voxel
template<unsigned short DIM, class C1, class C2, class FUNCTION>
CartesianGrid<DIM, C1> localConvertCartesian(const CartesianGrid<DIM, C2>& grid0, FUNCTION rule);
} // namespace auxi

} // namespace convertGrid
} // namespace vox
} // namespace merope

#include "../Grid/ConvertGrix.ixx"

#endif /* GRID_CONVERTGRIX_HXX_ */
