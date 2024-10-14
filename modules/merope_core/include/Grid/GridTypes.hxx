//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

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
using GridPhaseFrac = CartesianGrid<DIM, composite::Iso<PhaseType>>;

//! grid containing a single field
template<unsigned short DIM>
using GridField = CartesianGrid<DIM, double>;

//! grid containing phases
template<unsigned short DIM>
using GridPhase = CartesianGrid<DIM, PhaseType>;

}  // namespace vox
}  // namespace merope


