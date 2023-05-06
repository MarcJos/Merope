//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MEROPE_CORE_SRC_VOXELLATION_VOXRECURSTRUCTURE_IXX_
#define MEROPE_CORE_SRC_VOXELLATION_VOXRECURSTRUCTURE_IXX_


#include "VoxSimpleGauss.hxx"
#include "VoxSimpleMultiInclusions.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM, class VOXEL_TYPE, class BasicStruct>
inline CartesianGrid<DIM, VOXEL_TYPE> transformIntoGrid(
    const BasicStruct& basicStruct,
    const vox::PreSubGrid<DIM>& gridParameters) {
    //
    constexpr bool isMultiInclusions = std::is_same<BasicStruct, MultiInclusions<DIM>>::value;
    constexpr bool isField = std::is_same<BasicStruct, CartesianField<DIM>>::value;
    static_assert(isMultiInclusions or isField);
    //
    if constexpr (isMultiInclusions) {
        return VoxSimpleMultiInclusions<DIM, GetVoxelRule<VOXEL_TYPE>()>(basicStruct, gridParameters)();
    }
    if constexpr (isField) {
        return vox::VoxSimpleGauss<DIM>(basicStruct, gridParameters)();
    }
}

template <unsigned short DIM, class VOXEL_TYPE, class BasicStruct, class BasicType>
inline void VoxStructure<DIM, VOXEL_TYPE, BasicStruct, BasicType>::build() {
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

} //namespace vox
} // namespace merope


#endif /* MEROPE_CORE_SRC_VOXELLATION_VOXRECURSTRUCTURE_IXX_ */
