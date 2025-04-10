//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../GenericTools/Loops.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"

#include "../../../AlgoPacking/include/MultiDArrayObject.hxx"

#include "../SingleVoxel/SingleVoxel_Headers.hxx"
#include "../Grid/PreGrid.hxx"

namespace merope {
namespace vox {

template<unsigned short DIM, class VOXEL_TYPE>
//! general class for cartesian grid of voxels containing an information of a certain type
//! can use subgrid (price is perf)
class CartesianGrid : public MultiDArrayObject<DIM, VOXEL_TYPE, SubArrayDimensions<DIM>> {
public:
    //! constructor
    CartesianGrid(GridParameters<DIM> gridParameters_) :
        MultiDArrayObject<DIM, VOXEL_TYPE, SubArrayDimensions<DIM>>(
            SubArrayDimensions<DIM>(gridParameters_.getNbNodes(), gridParameters_.getNMin(), gridParameters_.getNMax()),
            VOXEL_TYPE{}
        ),
        gridParameters(gridParameters_) {}
    //! constructor
    CartesianGrid(GridParameters<DIM> gridParameters_, VOXEL_TYPE voxel) :
        MultiDArrayObject<DIM, VOXEL_TYPE, SubArrayDimensions<DIM>>(
            SubArrayDimensions<DIM>(gridParameters_.getNbNodes(), gridParameters_.getNMin(), gridParameters_.getNMax()),
            voxel),
        gridParameters(gridParameters_) {}
    //! getter
    const GridParameters<DIM>& getGridParameters() const { return gridParameters; }
    //! getter
    const array<size_t, DIM>& getNbNodeSubGrid() const { return this->getGridParameters().getNbNodeSubGrid(); }
    //! getter
    const array<size_t, DIM>& getNbNodeBigGrid() const { return this->getGridParameters().getNbNodeBigGrid(); }
    //! getter
    const array<size_t, DIM>& getNMin() const { return this->getGridParameters().getNMin(); }
    //! getter
    const array<size_t, DIM>& getNMax() const { return this->getGridParameters().getNMax(); }
    //! getter
    const Point<DIM>& getDx() const { return this->getGridParameters().getDx(); }
    //! getter
    const Point<DIM>& getL() const { return this->getGridParameters().getL(); }
    //! getter
    template<class INT_TYPE>
    Point<DIM> getCenterVoxel(const array<INT_TYPE, DIM>& ijk) const;

private:
    //! grid parameters. Never modify by setter!
    GridParameters<DIM> gridParameters;
};

} /* namespace vox */
}  // namespace merope

#include "../Grid/CartesianGrid.ixx"


