//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../SingleVoxel/GridTypesBase.hxx"


namespace merope {
namespace vox {
namespace composite {

//! @brief :  change only the phase of a composite
//! @param conversion : function turning a phase into another one (with possibly a different type)
template<unsigned short DIM, class PHASE_OUT, class COMPOSITE,
    class CONVERSION, typename = enable_if_t<is_composite<COMPOSITE>>>
Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> transform_phase(const COMPOSITE& local_data,
    const CONVERSION& conversion);

//! @brief : apply a texture on a composite
//! @param texturer : is a function depending on the phase + spatial position
template<unsigned short DIM, class COMPOSITE, class TEXTURER,
    class PHASE_OUT = double, typename = enable_if_t<is_composite<COMPOSITE>>>
Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> apply_texture_loc(const COMPOSITE& local_data,
    const TEXTURER& texturer, const Point<DIM>& x);


//! @brief : convert naturally from voxel types to others
template<unsigned short DIM, class COMPOSITE_OUT, class COMPOSITE_IN>
COMPOSITE_OUT convert_to(const COMPOSITE_IN& composite_in);

//! @return a function inserting inside a given voxel a pure phase (=assumed to cover the whole voxel)
//! @param purePhase : purePhase to be inserted
template<unsigned short DIM, class VOXEL_POLICY, class PHASE_TYPE>
void pure_phase_inserter(OutputFormat<VOXEL_POLICY::voxelRule, DIM, PHASE_TYPE>&, PHASE_TYPE purePhase, const VOXEL_POLICY& voxelPolicy);

//! @brief fill the voxel_to_include into the voxel_to_be_filled, according to the voxelPolicy
//! @param voxel_to_be_filled 
//! @param voxel_to_include 
//! @param voxelPolicy : @see vox::VoxelPolicy
template<unsigned short DIM, class VOXEL_TYPE, class VOXEL_POLICY>
void voxel_filler(VOXEL_TYPE& voxel_to_be_filled, const VOXEL_TYPE& voxel_to_include, const VOXEL_POLICY& voxelPolicy);

//! @brief : post-process the voxel to put it into a correct final state. (For example : the sum of volume fraction = 1)
template<unsigned short DIM, class VOXEL_TYPE, class VOXEL_POLICY>
void postProcess(VOXEL_TYPE& voxel, const VOXEL_POLICY& voxelPolicy, const auto& matrixPhaseHolder);

//! @brief : define a default voxel
template<class VOXEL_TYPE, class VOXEL_POLICY>
VOXEL_TYPE make_default_voxel(const VOXEL_POLICY& voxelPolicy, const auto& matrixPhaseHolder);


namespace auxi {
//! @brief : class for voxel conversion. Useful for partial specialization.
template<unsigned short DIM, class COMPOSITE_OUT, class COMPOSITE_IN>
struct Helper_convert_to {
    //!
    static COMPOSITE_OUT eval(const COMPOSITE_IN& composite_in);
private:
    //! @brief 
    template<class A, class B>
    static constexpr bool areTypes = is_same_v<COMPOSITE_OUT, A> and is_same_v<COMPOSITE_IN, B>;
};

//! @return the voxel written in vtk format
template<unsigned short DIM, class COMPOSITE>
to_stl_format<COMPOSITE> convert_to_stl_format(const COMPOSITE& composite_t);

}  // namespace  auxi

}  // namespace  composite
}  // namespace vox
}  // namespace merope

#include "ManipOn_GridTypeBases.ixx"


