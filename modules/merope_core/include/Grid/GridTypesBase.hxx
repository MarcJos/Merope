//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../Grid/CompositeVoxelTypes.hxx"
#include "../Grid/ListPhaseFrac.hxx"

#include "../MeropeNamespace.hxx"
#include "../Grid/VoxelRule.hxx"


namespace merope {
namespace vox {

//! representation of a single volume fraction of a phase
using SinglePhaseFrac = PhaseFrac<PhaseType>;

namespace composite {
//! @brief Class for defining 'pure' composite voxels
template<class TYPE_PHASE>
using Pure = TYPE_PHASE;

//! \briefClass for defining composite voxels
//! isotropic mixture of phases with different volume fractions
template<class TYPE_PHASE>
using Iso = gridAuxi::ListPhaseFrac<PhaseFrac<TYPE_PHASE>>;

//! \briefClass for defining composite voxels
//! mixture of phases with different volume fractions
//! with a direction of lamination
template<unsigned short DIM, class TYPE_PHASE>
using AnIso = gridAuxi::ListPhaseFrac<PhaseFracNormal<DIM, TYPE_PHASE>>;


//! @brief  whether is a vox::composite::Pure
template<class COMPOSITE>
constexpr bool is_Pure = is_arithmetic_v<COMPOSITE>;
//!  \brief whether is a list of PhaseFrac
template<class COMPOSITE>
constexpr bool is_Iso = false;
template<class PHASE_TYPE>
constexpr bool is_Iso<composite::Iso<PHASE_TYPE>> = true;
//! @brief  whether is a vox::composite::AnIso
template<class COMPOSITE>
constexpr bool is_AnIso = false;
template<unsigned short DIM, class PHASE_TYPE>
constexpr bool is_AnIso<vox::composite::AnIso<DIM, PHASE_TYPE>> = true;

template<class COMPOSITE>
constexpr bool is_composite = is_Iso<COMPOSITE> or is_AnIso<COMPOSITE> or is_Pure<COMPOSITE>;

template<class COMPOSITE, typename = std::enable_if_t<is_composite<COMPOSITE>>>
using Basic_Phase_Type = conditional_t<
    is_same_v<COMPOSITE, PhaseType>
    or is_same_v<COMPOSITE, Iso<PhaseType>>
    or is_same_v<COMPOSITE, AnIso<2, PhaseType>>
    or is_same_v<COMPOSITE, AnIso<3, PhaseType>>
    , PhaseType, double>;



template<unsigned short DIM, class COMPOSITE, class TYPE_PHASE, typename = std::enable_if_t<is_composite<COMPOSITE>>>
using Change_Type_Composite =
conditional_t<is_Pure<COMPOSITE>, Pure<TYPE_PHASE>,
    conditional_t<is_Iso<COMPOSITE>, Iso<TYPE_PHASE>,
    conditional_t<is_AnIso<COMPOSITE>, AnIso<DIM, TYPE_PHASE>,
    void*
    >>>;



template<vox::VoxelRule voxelRule, unsigned short DIM, class PHASE_TYPE>
using OutputFormat = typename std::conditional_t<voxelRule == vox::VoxelRule::Center, vox::composite::Pure<PHASE_TYPE>,
    std::conditional_t< voxelRule == vox::VoxelRule::Average, vox::composite::Iso<PHASE_TYPE>,
    vox::composite::AnIso<DIM, PHASE_TYPE>>>;

template<class TYPE_PHASE>
using stl_format_Pure = TYPE_PHASE;
template<class TYPE_PHASE>
using stl_format_Iso = vector<tuple<TYPE_PHASE, double>>;
template<unsigned short DIM, class TYPE_PHASE>
using stl_format_AnIso = tuple<stl_format_Iso<TYPE_PHASE>, Point<DIM>>;

template<class COMPOSITE, typename = std::enable_if_t<is_composite<COMPOSITE>>>
using to_stl_format =
conditional_t <
    is_Pure<COMPOSITE>, stl_format_Pure<Basic_Phase_Type<COMPOSITE>>,
    conditional_t<
    is_Iso<COMPOSITE>, stl_format_Iso<Basic_Phase_Type<COMPOSITE>>,
    conditional_t<
    is_same_v<COMPOSITE, AnIso<2, PhaseType>>, stl_format_AnIso<2, PhaseType>,
    conditional_t<
    is_same_v<COMPOSITE, AnIso<3, PhaseType>>, stl_format_AnIso<3, PhaseType>,
    conditional_t<
    is_same_v<COMPOSITE, AnIso<2, double>>, stl_format_AnIso<2, double>,
    conditional_t<
    is_same_v<COMPOSITE, AnIso<3, double>>, stl_format_AnIso<3, double>,
    void*
    >>>>>>;


//! @return the voxel written in vtk format
template<unsigned short DIM, class COMPOSITE>
to_stl_format<COMPOSITE> convert_to_stl_format(const COMPOSITE& composite_t);

}  // namespace  composite

template<class OUTPUT_FORMAT>
constexpr vox::VoxelRule GetVoxelRule() {
    if constexpr (std::is_same_v<OUTPUT_FORMAT, vox::PhaseType>) {
        return vox::VoxelRule::Center;
    } else if constexpr (std::is_same_v<OUTPUT_FORMAT, vox::composite::Iso<PhaseType>>) {
        return vox::VoxelRule::Average;
    } else {
        return vox::VoxelRule::Laminate;
    }
}

}  // namespace vox
}  // namespace merope

#include "GridTypesBase.ixx"


