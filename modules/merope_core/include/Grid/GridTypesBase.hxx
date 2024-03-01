//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef VOXELLATION_GRIDTYPEBASE_HXX_
#define VOXELLATION_GRIDTYPEBASE_HXX_

#include "../../../AlgoPacking/src/SphereContainer.hxx"
#include "../Grid/ListPhaseFrac.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

//! phase written in the .vtk files read by tmfft and amitex
using VTK_PHASE = unsigned short;
//! representation of a single volume fraction of a phase
using SinglePhaseFrac = auxi_SphereCollection::PhaseFrac<VTK_PHASE>;
//! representation of the content of a volume filled by different phases in different fractions
using VoxelPhaseFrac = gridAuxi::ListPhaseFrac<VTK_PHASE>;
//! representation of the content of a volume filled by different values in different fractions
using VoxelValueFrac = gridAuxi::ListPhaseFrac<double>;
//! whether is a list of PhaseFrac
template<class POTENTIAL_VOX_TYPE_FRAC>
constexpr bool is_Voxel_TYPE_Frac = (is_same<VoxelPhaseFrac, POTENTIAL_VOX_TYPE_FRAC>::value
    or is_same<VoxelValueFrac, POTENTIAL_VOX_TYPE_FRAC>::value);

} // namespace vox
} // namespace merope

#endif /* VOXELLATION_GRIDTYPEBASE_HXX_ */
