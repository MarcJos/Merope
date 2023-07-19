//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MEROPE_CORE_SRC_VOXELLATION_VOXRECURSTRUCTURE_HXX_
#define MEROPE_CORE_SRC_VOXELLATION_VOXRECURSTRUCTURE_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Grid/VoxelRule.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../Grid/PreGrid.hxx"
#include "../MesoStructure/Structure.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

//! Virtual class for building a (scalar) field
template<unsigned short DIM, class VOXEL_TYPE>
class VoxGrid {
    using GRID_TYPE = CartesianGrid<DIM, VOXEL_TYPE>;
public:
    VoxGrid(PreSubGrid<DIM> gridParameters_) :
        gridParameters{ gridParameters_ },
        voxGrid{ nullptr } {
        voxGrid = make_unique<GRID_TYPE>(GRID_TYPE(gridParameters));
    }
    //! destructor
    virtual ~VoxGrid() {}
    //! getter on the grid of volume fractions
    GRID_TYPE operator()() {
        build();
        return this->getVoxGrid();
    }
    //! getter
    const PreSubGrid<DIM>& getGridParameters() const { return gridParameters; }
protected:
    //! fills the vector phaseFracVol
    virtual void build() = 0;
    //! getter
    GRID_TYPE& getVoxGrid() { return *voxGrid; }
private:
    //! grid parameters
    PreSubGrid<DIM> gridParameters;
    //! vector describing on each voxel the volume fraction of a given initial phase
    unique_ptr<GRID_TYPE> voxGrid;
};

template<unsigned short DIM, class VOXEL_TYPE, class BasicStruct, class BasicType>
//! Class for building a phase field from a Structure geometry adapted to the recursive shape of RecursiveStructure
class VoxStructure : public VoxGrid<DIM, VOXEL_TYPE> {
    using SELF_TYPE = VoxStructure<DIM, VOXEL_TYPE, BasicStruct, BasicType>;
    using RecurStruct = RecursiveStructure<DIM, BasicStruct, BasicType>;
public:
    //! constructor
    VoxStructure(const RecurStruct& structure_, PreSubGrid<DIM> gridParameters_) :
        VoxGrid<DIM, VOXEL_TYPE>(gridParameters_), structure{ &structure_ } {}
protected:
    //! inner representation
    const RecurStruct* structure;
private:
    //! fills the vector phaseFracVol
    void build() override;
};

//! turns a basic structure (=Field, MultiInclusions) into a CartesianGrid of prescribed descriptions
//! \warning : for each basic structure, there is a natural way of transforming them
//! \param structure : input to voxellize
//! \param gridParameters : parameters of the grid
template<unsigned short DIM, class VOXEL_TYPE, class BasicStruct>
vox::CartesianGrid<DIM, VOXEL_TYPE> transformIntoGrid(const BasicStruct& basicStruct, const vox::PreSubGrid<DIM>& gridParameters);

} // namespace vox
} // namespace merope


#include "VoxRecurStructure.ixx"

#endif /* MEROPE_CORE_SRC_VOXELLATION_VOXRECURSTRUCTURE_HXX_ */
