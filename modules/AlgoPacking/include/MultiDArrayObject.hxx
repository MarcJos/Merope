//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "../../Geometry/include/GeomTypes.hxx"

#include "ArrayDimensions.hxx"

namespace merope {
namespace vox {

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
class MultiDArrayObject : public vector<VOXEL_TYPE> {
public:
    //! default constructor
    MultiDArrayObject() : vector<VOXEL_TYPE>{}, arrayDim(create_array<DIM, size_t>(0)) {}
    //! repetitive constructor
    MultiDArrayObject(ARRAY_DIMENSION arrayDim_, VOXEL_TYPE voxel)
        : vector<VOXEL_TYPE>(arrayDim_.getTotalNumberVoxels(), voxel), arrayDim(arrayDim_) {}

    //! fill all the voxels with the same data
    void fillAll(VOXEL_TYPE voxelData);
    //! fill a a slice of the grid with the data, i.e. with i and j fixed, and only k variable
    //! hypothesis : the projection of limits is indeed inside the subgrid
    //! \param ijk : point located on the line of the slice
    //! \param limits : extremities of the slice (in the last direction)
    //! \param data : data to be filled into the slice (always the same)
    template<class INDEX_TYPE, class LOCAL_FUNCTION>
    void apply_voxelSlice(array<INDEX_TYPE, DIM> ijk, array<INDEX_TYPE, 2> limits, const LOCAL_FUNCTION& local_function);
    //! using the classical operator []
    using vector<VOXEL_TYPE>::operator[];
    //! using the operator with entire indices (ijk is assumed to be within the grid). Beware, expensive!
    VOXEL_TYPE& operator[](const array<size_t, DIM>& ijk);
    //! using the operator with entire indices (ijk is assumed to be within the grid). Beware, expensive!
    const VOXEL_TYPE& operator[](const array<size_t, DIM>& ijk) const;
    //! using the operator with periodic indices. Beware, expensive!
    VOXEL_TYPE& operator[](const array<long, DIM>& ijk);
    //! using the operator with periodic indices. Beware, expensive!
    const VOXEL_TYPE& operator[](const array<long, DIM>& ijk) const;
    //! getter
    const ARRAY_DIMENSION& getArrayDimensions() const { return arrayDim; }

private:
    //! tests wheter an array is a subarray
    static constexpr bool isSubArray = is_same_v<ARRAY_DIMENSION, vox::SubArrayDimensions<DIM>>;
    //! return the linear index on the subGrid, assuming that ijk is is inside the grid
    //! \param ijk : discrete point to be projected onto the periodic grid
    size_t get_linear_index(const array<size_t, DIM>& ijk) const;
    //! return the linear index on the subGrid, assuming that ijk is is inside the grid
    //! \param ijk : discrete point to be projected onto the periodic grid
    size_t get_linear_index_periodic(const array<long, DIM>& ijk) const;
    //! in case of a subgrid, substracts the nMin
    array<size_t, DIM> renormalize_index(array<size_t, DIM>ijk) const;
    //! fill a voxel with the data (fast for available linearIndex)
    template<class LOCAL_FUNCTION>
    void apply_voxelLinear(size_t linearIndex, const LOCAL_FUNCTION& local_function) { local_function((*this)[linearIndex]); }
    //! array dimensions
    ARRAY_DIMENSION arrayDim;
};

namespace auxi {
//! \return the projection of the periodic index onto the discrete grid
//! \param ijk : discrete point to be projected onto the periodic grid
//! \param nbNodes : dimensions of the discrete periodic grid
template<unsigned short DIM, class C>
array<size_t, DIM> project_index_periodic(const array<long, DIM>& ijk, const C& nbNodes);
//! \return the linear index on the grid of nbNodes, assuming that ijk is not necessarily on the grid
//! \param ijk : discrete point to be projected onto the periodic grid
//! \param nbNodes : dimensions of the discrete periodic grid
template<unsigned short DIM, class C>
size_t get_linear_index_periodic(const array<long, DIM>& ijk, const C& nbNodes);
//! \return the linear index on the grid of nbNodes, assuming that ijk is is inside the grid
//! \param ijk : discrete point to be projected onto the periodic grid
//! \param nbNodes : dimensions of the discrete periodic grid
template<unsigned short DIM, class C>
size_t get_linear_index(const array<size_t, DIM>& ijk, const C& nbNodes);
//! \return : a list of cuboid on the square defined by nbNodes corresponding to the projection of a periodic cuboid
//! \param : cuboidPerLimits = {{imin, imax}, {jmin, jmax}, {kmin, kmax}} for ther periodic cuboid [imin, imax) x [jmin, jmax) x [kmin, kmax)
//! \param nbNodes : periodicity lengths (in each direction)
template<unsigned short DIM>
vector<array<array<size_t, 2>, DIM>> getCuboidList(const array<array<long, 2>, DIM>& cuboidPerLimits, const array<size_t, DIM>& nbNodes);
//! projet a segment [i, j) on the periodic torus of dimension nbNode onto the real line,
//! \return : the corresponding list of segments in [0,nbNode)
//! \param perLimits : {i, j} limits of the periodic segment
//! \param nbNode : periodicity length
vector<array<size_t, 2>> getSegmentsFromPeriodic(const array<long, 2>& perLimits, size_t nbNode);
//! loop for accessing rapidly the linear indices corresponding to triplets ijk inside the cuboid described by limits
//! for each of those indices linearIndex, apply f(linearIndex)
//! usefull for fast writing in cuboids
template<unsigned short DIM, class FUNCTION>
void loopOnIndex(const array<array<size_t, 2>, DIM>& limits, const array<size_t, DIM>& nbNodes, const FUNCTION& f);
//! fixme : slice in the last direction
template<unsigned short DIM, class FUNCTION>
void loopOnIndexSlice(array<size_t, DIM> ijk, const array<size_t, 2>& limits, const array<size_t, DIM>& nbNodes, const FUNCTION& f);
}  // namespace auxi

}  // namespace vox

}  // namespace merope

#include "MultiDArrayObject.ixx"


