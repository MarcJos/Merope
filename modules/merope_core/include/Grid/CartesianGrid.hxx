//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_CARTESIANGRID_HXX_
#define GRID_CARTESIANGRID_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"
#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../../../AlgoPacking/src/Loops.hxx"
#include "../../../AlgoPacking/src/MultiDArrayObject.hxx"

#include "../MeropeNamespace.hxx"
#include "../Grid/ListPhaseFrac.hxx"
#include "../Grid/PreGrid.hxx"

namespace merope {
namespace vox {

template<unsigned short DIM, class VOXEL_TYPE>
//! general class for cartesian grid of voxels containing an information of a certain type
//! can use subgrid (price is perf)
class CartesianGrid : public MultiDArrayObject<DIM, VOXEL_TYPE, SubArrayDimensions<DIM>> {
public:
    //! constructor
    CartesianGrid(PreSubGrid<DIM> gridParameters_) :
        MultiDArrayObject<DIM, VOXEL_TYPE, SubArrayDimensions<DIM>>(
            SubArrayDimensions<DIM>(gridParameters_.getNbNodes(), gridParameters_.getNMin(), gridParameters_.getNMax()),
            VOXEL_TYPE{}
        ),
        gridParameters(gridParameters_) {}
    //! getter
    const PreSubGrid<DIM>& getGridParameters() const { return gridParameters; }
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

private:
    //! grid parameters. Never modify by setter!
    PreSubGrid<DIM> gridParameters;
};

template<unsigned short DIM, class VOXEL_TYPE>
//! general class for cartesian grid of voxels containing an information of a certain type
//! cannot use subgrid here
class EfficientCartesianGrid : public MultiDArrayObject<DIM, VOXEL_TYPE, ArrayDimensions<DIM>>,
    public With_dx<DIM> {
public:
    //! @brief 
    //! @param gridParameters_ : describes dx and nb voxels per edge
    EfficientCartesianGrid(PreGrid<DIM> gridParameters_) :
        MultiDArrayObject<DIM, VOXEL_TYPE, ArrayDimensions<DIM>>(
            ArrayDimensions<DIM>(gridParameters_.getNbNodes()), VOXEL_TYPE{}
        ),
        With_dx<DIM>(gridParameters_.getDx()),
        L(gridParameters_.getL()) {}
    //! getter
    const Point<DIM>& getL()const { return L; }

private:
    //! L
    Point<DIM> L;
};


} /* namespace vox */
} // namespace merope

#include "../Grid/CartesianGrid.ixx"

#endif /* GRID_CARTESIANGRID_HXX_ */
