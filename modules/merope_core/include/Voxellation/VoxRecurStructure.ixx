//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "VoxSimpleGauss.hxx"
#include "VoxSimpleMultiInclusions.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"



#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {


template<unsigned short DIM, class VOXEL_TYPE, class STRUCTURE>
vox::CartesianGrid<DIM, VOXEL_TYPE> transformIntoGrid(const  STRUCTURE& structure, const vox::GridParameters<DIM>& gridParameters) {
    static_assert(is_same_v<STRUCTURE, MultiInclusions<DIM>> or
        is_same_v<STRUCTURE, CartesianField<DIM>> or
        is_same_v<STRUCTURE, Structure<DIM>> or
        is_same_v<STRUCTURE, FieldStructure<DIM>>);
    //
    if constexpr (is_same_v<STRUCTURE, MultiInclusions<DIM>>) {
        return VoxSimpleMultiInclusions<DIM, GetVoxelRule<VOXEL_TYPE>()>(structure, gridParameters)();
    }
    // 
    else if constexpr (is_same_v<STRUCTURE, CartesianField<DIM>>) {
        auto field = vox::VoxSimpleGauss<DIM>(structure, gridParameters)();
        if constexpr (std::is_same_v<VOXEL_TYPE, double>) {
            return field;
        } else if constexpr (std::is_same_v<VOXEL_TYPE, vox::composite::Iso<double>>) {
            return vox::convertGrid::fromPureToIso(field);
        } else if constexpr (std::is_same_v<VOXEL_TYPE, vox::composite::AnIso<DIM, PhaseType>>) {
            return vox::convertGrid::fromPureToAnIso(field);
        } else {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Unexpected");
        }
    }
    //
    else if constexpr (is_same_v<STRUCTURE, Structure<DIM>>) {
        return auxi::VoxStructure<DIM, VOXEL_TYPE, MultiInclusions<DIM>, PhaseType>(structure, gridParameters)();
    }
    //
    else if constexpr (is_same_v<STRUCTURE, FieldStructure<DIM>>) {
        return auxi::VoxStructure<DIM, VOXEL_TYPE, CartesianField<DIM>, double>(structure, gridParameters)();
    }
}

template <unsigned short DIM, class VOXEL_TYPE, class BasicStruct, class BasicType>
void auxi::VoxStructure<DIM, VOXEL_TYPE, BasicStruct, BasicType>::build() {
    switch (structure->getTypeOf()) {
    case(auxiMicroStructure::TypeOfCombination::Simple):
    {
        this->getVoxGrid() = vox::transformIntoGrid<DIM, VOXEL_TYPE>(structure->getBasicStructure(), this->getGridParameters());
        break;
    }
    case (auxiMicroStructure::TypeOfCombination::Combination2):
    {
        auto struct1 = structure->getRecurStructure1();
        auto struct2 = structure->getRecurStructure2();
        const auto& gridFracVol1 = SELF_TYPE(struct1, this->getGridParameters())();
        const auto& gridFracVol2 = SELF_TYPE(struct2, this->getGridParameters())();
        this->getVoxGrid() = gridAuxi::combineGridFunction<DIM, VOXEL_TYPE>(gridFracVol1, gridFracVol2, structure->getTransformFunction());
        break;
    }
    case (auxiMicroStructure::TypeOfCombination::Mask):
    {
        auto struct1 = structure->getRecurStructure1();
        auto struct2 = structure->getRecurStructure2();
        auto mask = structure->getMask();
        const auto& gridFracVol1 = SELF_TYPE(struct1, this->getGridParameters())();
        const auto& gridFracVol2 = SELF_TYPE(struct2, this->getGridParameters())();
        const auto& gridFracVolMask = SELF_TYPE(mask, this->getGridParameters())();
        this->getVoxGrid() = gridAuxi::combineGridMask<DIM, VOXEL_TYPE>(gridFracVol1, gridFracVol2, gridFracVolMask);
        break;
    }
    default:
        cerr << __PRETTY_FUNCTION__ << endl;
        throw invalid_argument("Not implemented yet!");
    }
}

}  // namespace vox
}  // namespace merope



