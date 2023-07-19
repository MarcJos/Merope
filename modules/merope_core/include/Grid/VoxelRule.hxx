//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_VOXELRULE_HXX_
#define GRID_VOXELRULE_HXX_


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

//! parameters for computing the voxellation
enum class VoxelRule {
    Average,    // compute the percentage of each phase inside the voxel
    Center      // evaluate the phase at the very center of the voxel
};

class WithVoxelRule {
    //! class from which other inherit
public:
    WithVoxelRule(VoxelRule voxelRule_) : voxelRule{ voxelRule_ } {}
    //! set the rule for choosing which phase is in the voxel
    void setVoxelRule(VoxelRule voxelRule_) { voxelRule = voxelRule_; }
protected:
    //! determines which rule to fill the voxel
    VoxelRule voxelRule;
};

template<vox::VoxelRule voxelRule>
using OutputFormat = typename std::conditional<voxelRule == vox::VoxelRule::Average, vox::VoxelPhaseFrac, vox::VTK_PHASE>::type;

template<class OUTPUT_FORMAT>
constexpr vox::VoxelRule GetVoxelRule() {
    if constexpr (std::is_same<OUTPUT_FORMAT, vox::VoxelPhaseFrac>::value) {
        return vox::VoxelRule::Average;
    } else {
        return vox::VoxelRule::Center;
    }
}

} //namespace vox
} // namespace merope

#endif /* GRID_VOXELRULE_HXX_ */
