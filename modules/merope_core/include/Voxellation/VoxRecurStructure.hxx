//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../SingleVoxel/SingleVoxel_Headers.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../Grid/PreGrid.hxx"
#include "../MesoStructure/Structure.hxx"

#include "VoxGrid.hxx"


namespace merope {
namespace vox {

//! turns a basic structure (=Field, MultiInclusions) into a CartesianGrid of prescribed descriptions
//! \warning : for each basic structure, there is a natural way of transforming them
//! \param structure : input to voxellize
//! \param gridParameters : parameters of the grid
template<unsigned short DIM, class VOXEL_POLICY, class STRUCTURE>
vox::CartesianGrid<DIM, composite::OutputFormat<VOXEL_POLICY::voxelRule, DIM, Phase_Type_From_Structure_Type<STRUCTURE>>>
transformIntoGrid(const STRUCTURE& structure, const vox::GridParameters<DIM>& gridParameters,
    VOXEL_POLICY voxelPolicy);

namespace auxi {

template<unsigned short DIM, class VOXEL_POLICY, class BasicStruct, class BasicType>
//! Class for building a phase field from a Structure geometry adapted to the recursive shape of RecursiveStructure
class VoxStructure : public VoxGrid<DIM, composite::OutputFormat<VOXEL_POLICY::voxelRule, DIM, BasicType>> {
    using VOXEL_TYPE = composite::OutputFormat<VOXEL_POLICY::voxelRule, DIM, BasicType>;
    using SELF_TYPE = VoxStructure<DIM, VOXEL_POLICY, BasicStruct, BasicType>;
    using RecurStruct = RecursiveStructure<DIM, BasicStruct, BasicType>;
public:
    //! constructor
    VoxStructure(const RecurStruct& structure_, GridParameters<DIM> gridParameters_,
        VOXEL_POLICY voxelPolicy_) :
        VoxGrid<DIM, VOXEL_TYPE>(gridParameters_), structure{ &structure_ }, voxelPolicy{ voxelPolicy_ }{}

    const VOXEL_POLICY& getVoxelPolicy() const { return this->voxelPolicy; }
protected:
    //! inner representation
    const RecurStruct* structure;
    VOXEL_POLICY voxelPolicy;
private:
    //! fills the vector phaseFracVol
    void build() override;
};


}  // namespace  auxi

}  // namespace vox
}  // namespace merope


#include "VoxRecurStructure.ixx"


