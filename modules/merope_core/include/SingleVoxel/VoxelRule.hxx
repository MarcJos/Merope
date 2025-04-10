//! Copyright : see license.txt
//!
//! \brief Parameters describing how a voxellation should be done, namely, 
//! how the inner content of a single  voxels is computed
//! and how it it represented.
//!

#pragma once

#include "../../../GenericMerope/StdHeaders.hxx"

namespace merope {
namespace vox {

//! parameters for computing the voxellation
enum class VoxelRule {
    Center,         //!< @brief : evaluate the phase at the very center of the voxel
    Average,        //!< @brief : compute the percentage of each phase inside the voxel
    Laminate,       //!< @brief : compute the percentage of each phase inside the voxel + normal
    PolyGeom       //!< @brief : retain an approximate geometry inside the voxel, made of planes
};

//! @brief : class for determining which voxel policy should be applied when computing a voxel
//! @tparam VOXEL_RULE : inner content of a voxel
//! @tparam Assume_no_Intersection : assume that the geometrical shapes are not intersecting each others (eg non-intersecting spheres to be voxelled)
//! @tparam PHASE_TYPE : phase type representing the content of a voxel
template<VoxelRule VOXEL_RULE, bool Assume_no_Intersection, class PHASE_TYPE = int>
class VoxelPolicy;

//! @brief default implementation of voxel policy Assume no Intersection
template<VoxelRule VOXEL_RULE, class PHASE_TYPE>
class VoxelPolicy<VOXEL_RULE, true, PHASE_TYPE> {
public:
    static constexpr bool Assume_no_Intersection = true;
    static constexpr VoxelRule voxelRule = VOXEL_RULE;
};

//! @brief More complex of voxel policy allows for intersections
//! then rules should be chosen for dealing with them
template<VoxelRule VOXEL_RULE, class PHASE_TYPE>
class VoxelPolicy<VOXEL_RULE, false, PHASE_TYPE> {
public:
    static constexpr bool Assume_no_Intersection = false;
    static constexpr VoxelRule voxelRule = VOXEL_RULE;

    //! @brief : constructor
    //! the resulting phase when intersecting shapes of phases i1 and i2 is intersectionRule_(i1, i2)
    //! @param intersectionRule_ : rule for intersecting phases
    //! @warning : all the possible intersections shall be defined!
    VoxelPolicy(std::function<PHASE_TYPE(PHASE_TYPE, PHASE_TYPE)> intersectionRule_) :
        intersectionRule(std::move(intersectionRule_)) {}
    //! @brief : constructor from a map
    //! the resulting phase when intersecting shapes of phases i1 and i2 is intersectionRule_as_map[(i1, i2)]
    //! @param intersectionRule_as_map: map describing how the intersections are computed
    VoxelPolicy(std::map<tuple<PHASE_TYPE, PHASE_TYPE>, PHASE_TYPE> intersectionRule_as_map) :
        intersectionRule([intersectionRule_as_map](auto i1, auto i2) {
        try {
            return intersectionRule_as_map.at({ i1, i2 });
        }
        catch (const std::out_of_range& e) {
            Merope_assert(false, "The map for intersections cannot be evaluated on the tuple (" + to_string(i1) + "," + to_string(i2) + ")");
            return i1; // never seen. Just for avoiding message.
        }
            }) {}
    //! @brief : getter
    const std::function<PHASE_TYPE(PHASE_TYPE, PHASE_TYPE)>& getIntersectionRule() const { return intersectionRule; }

private:
    //! the resulting phase when intersecting shapes of phases i1 and i2 is intersectionRule(i1, i2)
    std::function<PHASE_TYPE(PHASE_TYPE, PHASE_TYPE)> intersectionRule;
};


}  // namespace vox
}  // namespace merope


