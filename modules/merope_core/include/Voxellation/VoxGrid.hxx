
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


namespace merope {
namespace vox {

//! Virtual class for building a (scalar) field
template<unsigned short DIM, class VOXEL_TYPE>
class VoxGrid {
    using GRID_TYPE = CartesianGrid<DIM, VOXEL_TYPE>;
public:
    VoxGrid(GridParameters<DIM> gridParameters_) :
        voxGrid{ make_unique<GRID_TYPE>(GRID_TYPE(gridParameters_)) } {}
    VoxGrid(GridParameters<DIM> gridParameters_, VOXEL_TYPE vox) :
        voxGrid{ make_unique<GRID_TYPE>(GRID_TYPE(gridParameters_, vox)) } {}
    //! destructor
    virtual ~VoxGrid() {}
    //! getter on the grid of volume fractions
    GRID_TYPE operator()() {
        build();
        return this->getVoxGrid();
    }
    //! getter
    const GridParameters<DIM>& getGridParameters() const { return voxGrid->getGridParameters(); }
protected:
    //! fills the vector phaseFracVol
    virtual void build() = 0;
    //! getter
    GRID_TYPE& getVoxGrid() { return *voxGrid; }
private:
    //! vector describing on each voxel the volume fraction of a given initial phase
    unique_ptr<GRID_TYPE> voxGrid;
};

}  // namespace vox
}  // namespace merope