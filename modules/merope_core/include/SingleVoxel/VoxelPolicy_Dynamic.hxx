//! Copyright : see license.txt
//!
//! \brief For simplicity, creates dynamically a VoxelPolicy

#pragma once

#include "../../../GenericMerope/StdHeaders.hxx"
#include "VoxelRule.hxx"

namespace merope {
namespace vox {

using VoxelPolicy_Dynamic = std::variant <
    VoxelPolicy<VoxelRule::Center, true, PhaseType>,
    VoxelPolicy<VoxelRule::Center, false, PhaseType>,
    VoxelPolicy<VoxelRule::Center, true, double>,
    VoxelPolicy<VoxelRule::Center, false, double>,
    VoxelPolicy<VoxelRule::Center, true, int>,
    VoxelPolicy<VoxelRule::Center, false, int>,

    VoxelPolicy<VoxelRule::Average, true, PhaseType>,
    VoxelPolicy<VoxelRule::Average, false, PhaseType>,
    VoxelPolicy<VoxelRule::Average, true, double>,
    VoxelPolicy<VoxelRule::Average, false, double>,
    VoxelPolicy<VoxelRule::Average, true, int>,
    VoxelPolicy<VoxelRule::Average, false, int>,

    VoxelPolicy<VoxelRule::Laminate, true, PhaseType>,
    VoxelPolicy<VoxelRule::Laminate, false, PhaseType>,
    VoxelPolicy<VoxelRule::Laminate, true, double>,
    VoxelPolicy<VoxelRule::Laminate, false, double>,
    VoxelPolicy<VoxelRule::Laminate, true, int>,
    VoxelPolicy<VoxelRule::Laminate, false, int>,

    VoxelPolicy<VoxelRule::PolyGeom, true, PhaseType>,
    VoxelPolicy<VoxelRule::PolyGeom, false, PhaseType>,
    VoxelPolicy<VoxelRule::PolyGeom, true, double>,
    VoxelPolicy<VoxelRule::PolyGeom, false, double>,
    VoxelPolicy<VoxelRule::PolyGeom, true, int>,
    VoxelPolicy<VoxelRule::PolyGeom, false, int>
> ;

template<bool Assume_no_Intersection, class PHASE_TYPE, class ...ARGS>
VoxelPolicy_Dynamic create_voxel_policy_dynamic(VoxelRule voxelRule, ARGS... args) {
    switch (voxelRule) {
    case (VoxelRule::Center):
        return VoxelPolicy_Dynamic(VoxelPolicy<VoxelRule::Center, Assume_no_Intersection, PHASE_TYPE>(args...));
        break;
    case (VoxelRule::Average):
        return VoxelPolicy_Dynamic(VoxelPolicy<VoxelRule::Average, Assume_no_Intersection, PHASE_TYPE>(args...));
        break;
    case (VoxelRule::Laminate):
        return VoxelPolicy_Dynamic(VoxelPolicy<VoxelRule::Laminate, Assume_no_Intersection, PHASE_TYPE>(args...));
        break;
    case (VoxelRule::PolyGeom):
        return VoxelPolicy_Dynamic(VoxelPolicy<VoxelRule::PolyGeom, Assume_no_Intersection, PHASE_TYPE>(args...));
        break;
    default:
        Merope_assert(false, "Unknown voxel rule");
        throw runtime_error("Useless, but avoiding false <<return errors>>");
    }
}

template<class FUNCTION>
inline void apply_function_to_voxel_policy(FUNCTION function, VoxelPolicy_Dynamic voxelPolicy_Dynamic) {
    switch (voxelPolicy_Dynamic.index()) {
    case 0: function(std::get<VoxelPolicy<VoxelRule::Center, true, PhaseType>>(voxelPolicy_Dynamic)); break;
    case 1: function(std::get<VoxelPolicy<VoxelRule::Center, false, PhaseType>>(voxelPolicy_Dynamic)); break;
    case 2: function(std::get<VoxelPolicy<VoxelRule::Center, true, double>>(voxelPolicy_Dynamic)); break;
    case 3: function(std::get<VoxelPolicy<VoxelRule::Center, false, double>>(voxelPolicy_Dynamic)); break;
    case 4: function(std::get<VoxelPolicy<VoxelRule::Center, true, int>>(voxelPolicy_Dynamic)); break;
    case 5: function(std::get<VoxelPolicy<VoxelRule::Center, false, int>>(voxelPolicy_Dynamic)); break;

    case 6: function(std::get<VoxelPolicy<VoxelRule::Average, true, PhaseType>>(voxelPolicy_Dynamic)); break;
    case 7: function(std::get<VoxelPolicy<VoxelRule::Average, false, PhaseType>>(voxelPolicy_Dynamic)); break;
    case 8: function(std::get<VoxelPolicy<VoxelRule::Average, true, double>>(voxelPolicy_Dynamic)); break;
    case 9: function(std::get<VoxelPolicy<VoxelRule::Average, false, double>>(voxelPolicy_Dynamic)); break;
    case 10: function(std::get<VoxelPolicy<VoxelRule::Average, true, int>>(voxelPolicy_Dynamic)); break;
    case 11: function(std::get<VoxelPolicy<VoxelRule::Average, false, int>>(voxelPolicy_Dynamic)); break;

    case 12: function(std::get<VoxelPolicy<VoxelRule::Laminate, true, PhaseType>>(voxelPolicy_Dynamic)); break;
    case 13: function(std::get<VoxelPolicy<VoxelRule::Laminate, false, PhaseType>>(voxelPolicy_Dynamic)); break;
    case 14: function(std::get<VoxelPolicy<VoxelRule::Laminate, true, double>>(voxelPolicy_Dynamic)); break;
    case 15: function(std::get<VoxelPolicy<VoxelRule::Laminate, false, double>>(voxelPolicy_Dynamic)); break;
    case 16: function(std::get<VoxelPolicy<VoxelRule::Laminate, true, int>>(voxelPolicy_Dynamic)); break;
    case 17: function(std::get<VoxelPolicy<VoxelRule::Laminate, false, int>>(voxelPolicy_Dynamic)); break;

    case 18: function(std::get<VoxelPolicy<VoxelRule::PolyGeom, true, PhaseType>>(voxelPolicy_Dynamic)); break;
    case 19: function(std::get<VoxelPolicy<VoxelRule::PolyGeom, false, PhaseType>>(voxelPolicy_Dynamic)); break;
    case 20: function(std::get<VoxelPolicy<VoxelRule::PolyGeom, true, double>>(voxelPolicy_Dynamic)); break;
    case 21: function(std::get<VoxelPolicy<VoxelRule::PolyGeom, false, double>>(voxelPolicy_Dynamic)); break;
    case 22: function(std::get<VoxelPolicy<VoxelRule::PolyGeom, true, int>>(voxelPolicy_Dynamic)); break;
    case 23: function(std::get<VoxelPolicy<VoxelRule::PolyGeom, false, int>>(voxelPolicy_Dynamic)); break;
    }
}

}  // namespace  vox
}  // namespace  merope
