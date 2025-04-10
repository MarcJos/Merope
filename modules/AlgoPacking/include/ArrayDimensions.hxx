//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../GenericMerope/StdHeaders.hxx"
#include "../../GenericTools/AuxiFunctions.hxx"

#include "../../Geometry/include/GeomTypes.hxx"
#include "../../Geometry/include/Point.hxx"

namespace merope {

namespace vox {

template<unsigned short DIM>
class ArrayDimensions {
public:
    //! constructor
    explicit ArrayDimensions(array<size_t, DIM> nbNodes_) : nbNodes(nbNodes_) {}
    //! setter
    void setNbNodes(array<size_t, DIM> nbNodes_) { nbNodes = nbNodes_; }
    //! getter
    const array<size_t, DIM>& getNbNodes() const { return nbNodes; }
    //! getter, size of the grid covering the torus
    const array<size_t, DIM>& getNbNodeBigGrid() const { return nbNodes; }
    //! getter, size of the subgrid
    const array<size_t, DIM>& getNbNodeSubGrid() const { return nbNodes; }
    //! getter
    size_t getTotalNumberVoxels() const { return auxi_function::productOf<size_t>(getNbNodes()); }
private:
    //! number of voxels in each direction
    array<size_t, DIM> nbNodes;
};

template<unsigned short DIM>
class SubArrayDimensions {
public:
    //! @brief  constructor
    //! @param nbNodes_
    explicit SubArrayDimensions(array<size_t, DIM> nbNodes_);
    /// @brief constructor
    /// @param nbNodes_
    /// @param nMin_
    /// @param nMax_
    SubArrayDimensions(array<size_t, DIM> nbNodes_, array<size_t, DIM> nMin_, array<size_t, DIM> nMax_);
    //! test whether everything is coherent
    bool testCoherent() const;
    //! getter
    const array<size_t, DIM>& getNMax() const { return nMax; }
    //! getter
    const array<size_t, DIM>& getNMin() const { return nMin; }
    //! getter
    size_t getTotalNumberVoxels() const { return auxi_function::productOf<size_t>(getNbNodeSubGrid()); }
    //! getter, size of the subgrid
    const array<size_t, DIM>& getNbNodeSubGrid() const { return nbNodeSubgrid; }
    //! getter, size of the grid covering the torus
    const array<size_t, DIM>& getNbNodeBigGrid() const { return nbNodes; }
    //! \return whether the grid covers the whole torus or not
    const bool& doesCoverTorus() const { return coverTorus; }
    //! set nbNodes, and default nMin and nMax
    void setNbNodes(array<size_t, DIM> nbNodes_);
    //! set nMin and nMax
    void setSubGridIndices(array<size_t, DIM> nMin_, array<size_t, DIM> nMax_);
    //! getter
    const array<size_t, DIM>& getNbNodes() const { return nbNodes; }

private:
    //! \return whether the grid covers the whole torus or not
    bool computeCoverTorus() const;
    //! \return the size of the subgrid
    void setNodeSubGrid();

private:
    //! number of voxels in each direction
    array<size_t, DIM> nbNodes;
    //! says whether the grid covers the whole torus or not
    bool coverTorus;
    // discrete coordinates of the lower corner
    array<size_t, DIM> nMin;
    // dicrete coorindates of the upper corner
    array<size_t, DIM> nMax;
    // nMax - nMin
    array<size_t, DIM> nbNodeSubgrid;
};

}  // namespace vox
}  // namespace merope

#include "ArrayDimensions.ixx"

