//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


namespace merope {
namespace vox {

template<unsigned short DIM, class VOXEL_POLICY>
VoxSimpleMultiInclusions<DIM, VOXEL_POLICY>::VoxSimpleMultiInclusions(
    const MultiInclusions<DIM>& multiI,
    GridParameters<DIM> gridParameters_,
    VOXEL_POLICY voxelPolicy_) :
    VoxGrid<DIM, vox::composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>>(gridParameters_, vox::composite::make_default_voxel<vox::composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>>(voxelPolicy_, MatrixPhaseHolder<PhaseType>(multiI))),
    MatrixPhaseHolder<PhaseType>(multiI),
    multiInclusions{ &multiI },
    halfDiagVoxel{ gridParameters_.getHalfDiagVoxel() },
    voxelPolicy(voxelPolicy_) {}

template<unsigned short DIM, class VOXEL_POLICY>
void VoxSimpleMultiInclusions<DIM, VOXEL_POLICY>::postProcessGrid() {
    const auto& matrixPhaseHolder_ = static_cast<MatrixPhaseHolder<PhaseType>>(*this);
    const auto& voxelPolicy_ = this->voxelPolicy;
    convertGrid::apply_inplace(this->getVoxGrid(),
        [&matrixPhaseHolder_, voxelPolicy_](auto& voxel) {
            composite::postProcess<DIM>(voxel, voxelPolicy_, matrixPhaseHolder_);
        }
    );
}

template<unsigned short DIM, class VOXEL_POLICY>
void VoxSimpleMultiInclusions<DIM, VOXEL_POLICY>::build() {
    multiInclusions->apply_on_all(
        [this](const auto& vecInclusions) {
            this->buildGridPhaseFrac_T_auxi(vecInclusions);
        }
    );
    this->postProcessGrid();
}

template<unsigned short DIM, class VOXEL_POLICY>
template<class C>
void vox::VoxSimpleMultiInclusions<DIM, VOXEL_POLICY>::buildGridPhaseFrac_T_auxi(const vector<C>& inclusions) {
    // would be nice to be parallel
    for (const auto& inclus : inclusions) {
        buildGridPhaseFrac_singleInclusion<C>(inclus);
    }
}

template<unsigned short DIM, class VOXEL_POLICY>
template<class INCLUSION_TYPE>
void vox::VoxSimpleMultiInclusions<DIM, VOXEL_POLICY>::buildGridPhaseFrac_singleInclusion(
    const INCLUSION_TYPE& inclusion) {
    auto torusGridLimits = vox::auxi::computeGridLimits<DIM>(inclusion.getCuboid(), this->getGridParameters());
    auto allGridLimits = vox::auxi::intersectGridLimits<DIM>(torusGridLimits, this->getGridParameters());
    for (const auto& gridLimits : allGridLimits) {
        vox::auxi::ConvexGrid<DIM> convexGrid{ gridLimits, this->getGridParameters().getDx() };
        convexGrid.template compute<VOXEL_POLICY::voxelRule>(inclusion);
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

template<unsigned short DIM, class VOXEL_POLICY>
template<class INCLUSION_TYPE>
void vox::VoxSimpleMultiInclusions<DIM, VOXEL_POLICY>::execute(vox::auxi::SliceInstruction<long> sliceInstrunction, DiscPoint<DIM> ijk,
    const INCLUSION_TYPE& inclusion) {
    if (sliceInstrunction.voxelType == vox::auxi::VoxelType::MonoPhase) {
        auto phase_to_insert = sliceInstrunction.phase;
        const auto& voxelPolicy_ = this->voxelPolicy;
        this->getVoxGrid().template apply_voxelSlice<long>(ijk, sliceInstrunction.limits,
            [phase_to_insert, &voxelPolicy_](auto& vox) {
                vox::composite::pure_phase_inserter<DIM>(vox, phase_to_insert, voxelPolicy_);
            });
    } else {  // (sliceInstrunction.voxelType == vox::auxi::VoxelType::Composite)
        for (ijk[DIM - 1] = sliceInstrunction.limits[0]; ijk[DIM - 1] < sliceInstrunction.limits[1]; ijk[DIM - 1]++) {
            computeInclusionVoxel(inclusion, ijk);
        }
    }
}

template<unsigned short DIM, class VOXEL_POLICY>
template<class C>
void VoxSimpleMultiInclusions<DIM, VOXEL_POLICY>::computeInclusionVoxel(const C& microInclusion, const DiscPoint<DIM>& indexVoxel) {
    auto& voxelData = this->getVoxGrid()[indexVoxel];
    vox::inclusionAndVoxel::fillVoxel<DIM, C, composite::OutputFormat<VOXEL_POLICY::voxelRule, DIM, PhaseType>>(
        microInclusion, indexVoxel, this->getGridParameters().getDx(), halfDiagVoxel, voxelData,
        voxelPolicy);
}

}  // namespace vox
}  // namespace merope
