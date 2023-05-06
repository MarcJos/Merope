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

template<class VOXEL_TYPE>
VOXEL_TYPE makeVoxelData(VTK_PHASE phase, double volumeFraction = 1.) {
    static_assert(is_same<VOXEL_TYPE, VTK_PHASE>::value or is_same<VOXEL_TYPE, VoxelPhaseFrac>::value);
    if constexpr (std::is_same<VOXEL_TYPE, VTK_PHASE>::value) {
        return phase;
    }
    else if constexpr (std::is_same<VOXEL_TYPE, VoxelPhaseFrac>::value) {
        return VoxelPhaseFrac({ SinglePhaseFrac(phase,volumeFraction) });
    }
}

} // namespace vox
} // namespace merope

#endif /* VOXELLATION_GRIDTYPEBASE_HXX_ */
