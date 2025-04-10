//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../MultiInclusions/MultiInclusions.hxx"
#include "VoxSimpleGauss.hxx"
#include "VoxSimpleMultiInclusions.hxx"


namespace merope {
namespace vox {


template<unsigned short DIM, class VOXEL_POLICY, class STRUCTURE>
vox::CartesianGrid<DIM, composite::OutputFormat<VOXEL_POLICY::voxelRule, DIM, Phase_Type_From_Structure_Type<STRUCTURE>>>
transformIntoGrid(const STRUCTURE& structure, const vox::GridParameters<DIM>& gridParameters,
    VOXEL_POLICY voxelPolicy) {
    using VOXEL_TYPE = composite::OutputFormat<VOXEL_POLICY::voxelRule, DIM, Phase_Type_From_Structure_Type<STRUCTURE>>;
    static_assert(is_same_v<STRUCTURE, MultiInclusions<DIM>> or
        is_same_v<STRUCTURE, CartesianField<DIM>> or
        is_same_v<STRUCTURE, Structure<DIM>> or
        is_same_v<STRUCTURE, FieldStructure<DIM>>);
    if constexpr (is_same_v<STRUCTURE, MultiInclusions<DIM>>) {
        return VoxSimpleMultiInclusions<DIM, VOXEL_POLICY>(structure, gridParameters, voxelPolicy)();
    }
    // 
    else if constexpr (is_same_v<STRUCTURE, CartesianField<DIM>>) {
        auto field = vox::VoxSimpleGauss<DIM>(structure, gridParameters)();
        if constexpr (std::is_same_v<VOXEL_TYPE, double>) {
            return field;
        } else if constexpr (std::is_same_v<VOXEL_TYPE, vox::composite::Iso<double>>
            or std::is_same_v<VOXEL_TYPE, vox::composite::AnIso<DIM, PhaseType>>) {
            return vox::convertGrid::convert_to<DIM, VOXEL_TYPE>(field);
        } else {
            Merope_assert(false, "Unexpected type!");
        }
    }
    //
    else if constexpr (is_same_v<STRUCTURE, Structure<DIM>>) {
        return auxi::VoxStructure<DIM, VOXEL_POLICY, MultiInclusions<DIM>, PhaseType>(structure, gridParameters, voxelPolicy)();
    }
    //
    else if constexpr (is_same_v<STRUCTURE, FieldStructure<DIM>>) {
        return auxi::VoxStructure<DIM, VOXEL_POLICY, CartesianField<DIM>, double>(structure, gridParameters, voxelPolicy)();
    }
}

template <unsigned short DIM, class VOXEL_TYPE, class BasicStruct, class BasicType>
void auxi::VoxStructure<DIM, VOXEL_TYPE, BasicStruct, BasicType>::build() {
    switch (structure->getTypeOf()) {
    case(auxiMicroStructure::TypeOfCombination::Simple):
    {
        this->getVoxGrid() = vox::transformIntoGrid<DIM>(structure->getBasicStructure(), this->getGridParameters(), this->getVoxelPolicy());
        break;
    }
    case (auxiMicroStructure::TypeOfCombination::Combination2):
    {
        auto struct1 = structure->template getRecurStructure<0>();
        auto struct2 = structure->template getRecurStructure<1>();
        const auto& gridFracVol1 = SELF_TYPE(struct1, this->getGridParameters(), this->getVoxelPolicy())();
        const auto& gridFracVol2 = SELF_TYPE(struct2, this->getGridParameters(), this->getVoxelPolicy())();
        this->getVoxGrid() = gridAuxi::combineGridFunction<DIM, VOXEL_TYPE>(gridFracVol1, gridFracVol2, structure->getTransformFunction());
        break;
    }
    case (auxiMicroStructure::TypeOfCombination::Mask):
    {
        auto struct1 = structure->template getRecurStructure<0>();
        auto struct2 = structure->template getRecurStructure<1>();
        auto mask = structure->getMask();
        const auto& gridFracVol1 = SELF_TYPE(struct1, this->getGridParameters(), this->getVoxelPolicy())();
        const auto& gridFracVol2 = SELF_TYPE(struct2, this->getGridParameters(), this->getVoxelPolicy())();
        const auto& gridFracVolMask = SELF_TYPE(mask, this->getGridParameters(), this->getVoxelPolicy())();
        this->getVoxGrid() = gridAuxi::combineGridMask<DIM, VOXEL_TYPE>(gridFracVol1, gridFracVol2, gridFracVolMask);
        break;
    }
    default:
        Merope_error_not_done();
    }
}

}  // namespace vox
}  // namespace merope
