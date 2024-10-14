//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Grid/VoxelRule.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../Grid/PreGrid.hxx"
#include "../MesoStructure/Structure.hxx"


#include "../MeropeNamespace.hxx"
#include "VoxGrid.hxx"


namespace merope {
namespace vox {

//! turns a basic structure (=Field, MultiInclusions) into a CartesianGrid of prescribed descriptions
//! \warning : for each basic structure, there is a natural way of transforming them
//! \param structure : input to voxellize
//! \param gridParameters : parameters of the grid
template<unsigned short DIM, class VOXEL_TYPE, class STRUCTURE>
vox::CartesianGrid<DIM, VOXEL_TYPE> transformIntoGrid(const STRUCTURE& structure, const vox::GridParameters<DIM>& gridParameters);

namespace auxi {

template<unsigned short DIM, class VOXEL_TYPE, class BasicStruct, class BasicType>
//! Class for building a phase field from a Structure geometry adapted to the recursive shape of RecursiveStructure
class VoxStructure : public VoxGrid<DIM, VOXEL_TYPE> {
    using SELF_TYPE = VoxStructure<DIM, VOXEL_TYPE, BasicStruct, BasicType>;
    using RecurStruct = RecursiveStructure<DIM, BasicStruct, BasicType>;
public:
    //! constructor
    VoxStructure(const RecurStruct& structure_, GridParameters<DIM> gridParameters_) :
        VoxGrid<DIM, VOXEL_TYPE>(gridParameters_), structure{ &structure_ } {}
protected:
    //! inner representation
    const RecurStruct* structure;
private:
    //! fills the vector phaseFracVol
    void build() override;
};


}  // namespace  auxi

}  // namespace vox
}  // namespace merope


#include "VoxRecurStructure.ixx"


