//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
inline VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::VoxSimpleMultiInclusions(const MultiInclusions<DIM>& multiI, GridParameters<DIM> gridParameters_) :
    VoxGrid<DIM, vox::composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>>(gridParameters_),
    multiInclusions{ &multiI }, matrixPresence{ multiI.is_there_matrix() }, matrixPhase{ multiI.getMatrixPhase() },
    halfDiagVoxel{ gridParameters_.getHalfDiagVoxel() }{
    if constexpr (VOXEL_RULE == VoxelRule::Center) {
        this->getVoxGrid().fillAll(matrixPhase);
    }
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
inline void VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::postProcessGrid() {
    if constexpr (VOXEL_RULE == VoxelRule::Average or VOXEL_RULE == VoxelRule::Laminate) {
        const auto matrixPhase_copy = matrixPhase;
        const auto is_there_matrix_copy = matrixPresence;
        convertGrid::apply_inplace(this->getVoxGrid(),
            [is_there_matrix_copy, matrixPhase_copy](auto& phfv) {
                assert(is_there_matrix_copy or phfv.size() != 0);
                phfv.merge();
                phfv.renormalize(is_there_matrix_copy, matrixPhase_copy);
            });
    }
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
inline void VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::build() {
    multiInclusions->apply_on_all(
        [this](const auto& vecInclusions) {
            this->buildGridPhaseFrac_T_auxi(vecInclusions);
        }
    );
    this->postProcessGrid();
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
template<class C>
inline void vox::VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::buildGridPhaseFrac_T_auxi(const vector<C>& inclusions) {
    // would be nice to be parallel
    for (const auto& inclus : inclusions) {
        buildGridPhaseFrac_singleInclusion<C>(inclus);
    }
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
template<class INCLUSION_TYPE>
inline void vox::VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::buildGridPhaseFrac_singleInclusion(
    const INCLUSION_TYPE& inclusion) {
    auto torusGridLimits = vox::auxi::computeGridLimits<DIM>(inclusion.getCuboid(), this->getGridParameters());
    auto allGridLimits = vox::auxi::intersectGridLimits<DIM>(torusGridLimits, this->getGridParameters());
    for (const auto& gridLimits : allGridLimits) {
        vox::auxi::ConvexGrid<DIM> convexGrid{ gridLimits, this->getGridParameters().getDx() };
        convexGrid.template compute<VOXEL_RULE>(inclusion);
#pragma omp parallel for // all ijk are independent
        for (size_t index = 0; index < convexGrid.getNbFirstIndices(); index++) {
            auto ijk = convexGrid.getIndexAllCoordinates(index);
            auto allSliceInstructions = convexGrid.getSliceInstructions(index);
            for (size_t i_0 = 0; i_0 < allSliceInstructions.size(); i_0++) {
                execute(allSliceInstructions[i_0], ijk, inclusion);
            }
        }
    }
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
template<class INCLUSION_TYPE>
inline void vox::VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::execute(vox::auxi::SliceInstruction<long> sliceInstrunction, DiscPoint<DIM> ijk,
    const INCLUSION_TYPE& inclusion) {
    if (sliceInstrunction.voxelType == vox::auxi::VoxelType::MonoPhase) {
        this->getVoxGrid().template fillSlice<long>(ijk, sliceInstrunction.limits,
            vox::composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>(sliceInstrunction.phase));
    } else {  // (sliceInstrunction.voxelType == vox::auxi::VoxelType::Composite)
        for (ijk[DIM - 1] = sliceInstrunction.limits[0]; ijk[DIM - 1] < sliceInstrunction.limits[1]; ijk[DIM - 1]++) {
            computeInclusionVoxel(inclusion, ijk);
        }
    }
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
template<class C>
void VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::computeInclusionVoxel(const C& microInclusion, const DiscPoint<DIM>& indexVoxel) {
    auto& voxelData = this->getVoxGrid()[indexVoxel];
    vox::inclusionAndVoxel::fillVoxel<DIM, C, composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>>(microInclusion, indexVoxel, this->getGridParameters().getDx(), halfDiagVoxel, voxelData);
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
void VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::fillVoxel(const DiscPoint<DIM>& indexVoxel, const vox::composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>& phases2Include) {
    auto& voxelData = this->getVoxGrid()[indexVoxel];
    vox::inclusionAndVoxel::fillVoxel<composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>>(voxelData, phases2Include);
}

}  // namespace vox
}  // namespace merope



