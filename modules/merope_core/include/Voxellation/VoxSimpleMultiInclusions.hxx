//! Copyright : see license.txt
//!
//! \briefDetails of implementation for the voxellisation
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Geometry/GeomTools.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../Grid/ConvertGrix.hxx"
#include "../Grid/ConvexGrid.hxx"
#include "../Grid/GridManipulations.hxx"
#include "../Grid/GridTypes.hxx"
#include "../Grid/InclusionAndVoxel.hxx"
#include "../Grid/PreGrid.hxx"
#include "../Grid/VoxelRule.hxx"
#include "../MultiInclusions/LaguerreTess.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"
#include "../MultiInclusions/SphereInclusions.hxx"
#include "../MesoStructure/Structure.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

//! Class for building a phase field from a MultiInclusions geometry
template<unsigned short DIM, vox::VoxelRule VOXEL_RULE>
class VoxSimpleMultiInclusions : public VoxGrid<DIM, vox::composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>> {
    // Preconditions
    static_assert(VOXEL_RULE == VoxelRule::Center or VOXEL_RULE == VoxelRule::Average or VOXEL_RULE == VoxelRule::Laminate);
    //
public:
    //! main constructor
    VoxSimpleMultiInclusions(const MultiInclusions<DIM>& multiI, GridParameters<DIM> gridParameters);

protected:
    //! inner inclusions ref
    const MultiInclusions<DIM>* multiInclusions;
    //! is there a matrix ?
    bool matrixPresence;
    //! phase of the matrix
    PhaseType matrixPhase;
    //! fills the vector phaseFracVol
    void build() override;

private:
    //! stores the value of the half-diagonal of a voxel
    double halfDiagVoxel;
    //! inspects whether the microInclusion intersects the voxel. If yes, fills it accordinng the VoxelRule
    //! \param indexVoxel : 3-D position of the voxel
    //! \param microInclusion : microInclusion
    //! \param dx : dimensions of the voxel
    template<class C>
    void computeInclusionVoxel(const C& microInclusion, const DiscPoint<DIM>& indexVoxel);
    //! simple method for filling a voxel with a given list of phases and volume fractions
    void fillVoxel(const DiscPoint<DIM>& indexVoxel, const vox::composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>& phases2Include);
    //! fills the vector phaseFracVol, considering a single class of inclusions
    template<class C>
    void buildGridPhaseFrac_T_auxi(const vector<C>& inclusions);
    //! modifies the vector phaseFracVol, considering a single inclusion
    //! considering a single inclusion, fills a cartesian grid by identifying the discrete boundary of the inclusion
    //! \param inclusion : the considered inclusion
    template<class C>
    void buildGridPhaseFrac_singleInclusion(const C& inclusions);
    //! guarantee that the grid has the desired properties (for phaseFracVol, each voxel is filled at 100%)
    void postProcessGrid();
    //! execute a single slice instruction
    template<class INCLUSION_TYPE>
    void execute(vox::auxi::SliceInstruction<long>, DiscPoint<DIM> ijk, const INCLUSION_TYPE& inclusion);
};

}  // namespace vox
}  // namespace merope


#include "VoxSimpleMultiInclusions.ixx"


