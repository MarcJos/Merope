//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_GRIDTYPES_HXX_
#define GRID_GRIDTYPES_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Grid/CartesianGrid.hxx"
#include "../Grid/GridTypesBase.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

////////////////////////////////////////////////////////////
///// All the types of grid
////////////////////////////////////////////////////////////

//! Cartesian grids
//! grid for phase field
template<unsigned short DIM>
using GridPhaseFrac = CartesianGrid<DIM, VoxelPhaseFrac>;

//! grid containing a single field
template<unsigned short DIM>
using GridField = CartesianGrid<DIM, double>;

//! grid containing phases
template<unsigned short DIM>
using GridPhase = CartesianGrid<DIM, VTK_PHASE>;

} // namespace vox
} // namespace merope

#endif /* GRID_GRIDTYPES_HXX_ */
