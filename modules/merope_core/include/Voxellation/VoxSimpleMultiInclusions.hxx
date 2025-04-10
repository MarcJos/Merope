//! Copyright : see license.txt
//!
//! \briefDetails of implementation for the voxellisation
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/GeomTools.hxx"

#include "../Grid/CartesianGrid.hxx"
#include "../Grid/ConvertGrix.hxx"
#include "../Grid/ConvexGrid.hxx"
#include "../Grid/GridManipulations.hxx"
#include "../Grid/GridTypes.hxx"
#include "../SingleVoxel/SingleVoxel_Headers.hxx"
#include "../Grid/PreGrid.hxx"
#include "../MultiInclusions/LaguerreTess.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"
#include "../MultiInclusions/SphereInclusions.hxx"
#include "../MesoStructure/Structure.hxx"
#include "../SingleVoxel/MatrixPhaseHolder.hxx"


namespace merope {
namespace vox {

//! Class for building a phase field from a MultiInclusions geometry
template<unsigned short DIM, class VOXEL_POLICY>
class VoxSimpleMultiInclusions : public VoxGrid<DIM, vox::composite::OutputFormat<VOXEL_POLICY::voxelRule, DIM, PhaseType>>, protected MatrixPhaseHolder<PhaseType> {
public:
    static constexpr vox::VoxelRule VOXEL_RULE = VOXEL_POLICY::voxelRule;
    //! main constructor
    VoxSimpleMultiInclusions(const MultiInclusions<DIM>& multiI, GridParameters<DIM> gridParameters,
        VOXEL_POLICY voxelPolicy);

protected:
    //! inner inclusions ref
    const MultiInclusions<DIM>* multiInclusions;
    //! fills the vector phaseFracVol
    void build() override;

private:
    //! inspects whether the microInclusion intersects the voxel. If yes, fills it accordinng the VoxelRule
    //! \param indexVoxel : 3-D position of the voxel
    //! \param microInclusion : microInclusion
    //! \param dx : dimensions of the voxel
    template<class C>
    void computeInclusionVoxel(const C& microInclusion, const DiscPoint<DIM>& indexVoxel);
    //! fills the vector phaseFracVol, considering a single class of inclusions
    template<class C>
    void buildGridPhaseFrac_T_auxi(const vector<C>& inclusions);
    //! modifies the vector phaseFracVol, considering a single inclusion
    //! considering a single inclusion, fills a cartesian grid by identifying the discrete boundary of the inclusion
    //! \param inclusion : the considered inclusion
    template<class C>
    void buildGridPhaseFrac_singleInclusion(const C& inclusions);
    //! @brief guarantee that the grid has the desired properties (for phaseFracVol, each voxel is filled at 100%)
    void postProcessGrid();
    //! execute a single slice instruction
    template<class INCLUSION_TYPE>
    void execute(vox::auxi::SliceInstruction<long>, DiscPoint<DIM> ijk, const INCLUSION_TYPE& inclusion);

private:
    //! stores the value of the half-diagonal of a voxel
    double halfDiagVoxel;
    //! special rule for filling voxels
    VOXEL_POLICY voxelPolicy;
};

}  // namespace vox
}  // namespace merope


#include "VoxSimpleMultiInclusions.ixx"


