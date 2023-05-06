//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef VOXELLATIONAUX_IXX_
#define VOXELLATIONAUX_IXX_


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
inline VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::VoxSimpleMultiInclusions(const MultiInclusions<DIM>& multiI, PreSubGrid<DIM> gridParameters_):
    VoxGrid<DIM, vox::OutputFormat<VOXEL_RULE>>(gridParameters_),
    multiInclusions{ &multiI }, matrixPhase{ multiI.getMatrixPhase() },
    halfDiagVoxel{ gridParameters_.getHalfDiagVoxel() }{
    if constexpr (VOXEL_RULE == VoxelRule::Center) {
        this->getVoxGrid().fillAll(matrixPhase);
    }
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
inline void VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::postProcessGrid() {
    if constexpr (VOXEL_RULE == VoxelRule::Average) {
        const auto matrixPhase_copy = matrixPhase;
        convertGrid::apply_inplace(this->getVoxGrid(), [matrixPhase_copy](auto& phfv) {
            phfv.merge();
            phfv.renormalize(matrixPhase_copy);
            return;
            });
    }
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
inline void VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::build() {
    buildGridPhaseFrac_T_auxi<smallShape::SphereInc<DIM>>(multiInclusions->getSphereInc());
    buildGridPhaseFrac_T_auxi<smallShape::ConvexPolyhedronInc<DIM>>(multiInclusions->getPolyhedrons());
    buildGridPhaseFrac_T_auxi<smallShape::EllipseInc<DIM>>(multiInclusions->getEllipseInc());
    buildGridPhaseFrac_T_auxi<smallShape::SpheroPolyhedronInc<DIM>>(multiInclusions->getSpheroPolyhedrons());
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
    INCLUSION_TYPE inclusion) {
    if constexpr (VOXEL_RULE == VoxelRule::Center) {
        if (sliceInstrunction.voxelType == vox::auxi::VoxelType::MonoPhase) {
            this->getVoxGrid().template fillSlice<long>(ijk, sliceInstrunction.limits, sliceInstrunction.phase);
        }
        else { // (sliceInstrunction.voxelType == vox::auxi::VoxelType::Composite)
            for (ijk[DIM - 1] = sliceInstrunction.limits[0]; ijk[DIM - 1] < sliceInstrunction.limits[1]; ijk[DIM - 1]++) {
                computeInclusionVoxel(inclusion, ijk);
            }

        }
    }
    if constexpr (VOXEL_RULE == VoxelRule::Average) {
        if (sliceInstrunction.voxelType == vox::auxi::VoxelType::MonoPhase) {
            this->getVoxGrid().template fillSlice<long>(ijk, sliceInstrunction.limits, VoxelPhaseFrac({ SinglePhaseFrac(sliceInstrunction.phase, 1.) }));
        }
        else { // (sliceInstrunction.voxelType == vox::auxi::VoxelType::Composite)
            for (ijk[DIM - 1] = sliceInstrunction.limits[0]; ijk[DIM - 1] < sliceInstrunction.limits[1]; ijk[DIM - 1]++) {
                computeInclusionVoxel(inclusion, ijk);
            }
        }
    }
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
template<class C>
void VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::computeInclusionVoxel(const C& microInclusion, const DiscPoint<DIM>& indexVoxel) {
    auto& voxelData = this->getVoxGrid()[indexVoxel];
    vox::inclusionAndVoxel::fillVoxel<DIM, C, OutputFormat<VOXEL_RULE>>(microInclusion, indexVoxel, this->getGridParameters().getDx(), halfDiagVoxel, voxelData);
}

template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
void VoxSimpleMultiInclusions<DIM, VOXEL_RULE>::fillVoxel(const DiscPoint<DIM>& indexVoxel, const vox::OutputFormat <VOXEL_RULE>& phases2Include) {
    auto& voxelData = this->getVoxGrid()[indexVoxel];
    vox::inclusionAndVoxel::fillVoxel<OutputFormat<VOXEL_RULE>>(voxelData, phases2Include);
}

} // namespace vox
} // namespace merope


#endif /* VOXELLATIONAUX_IXX_ */
