//! Copyright : see license.txt
//!
//! \brief

#pragma once


#include "Voxel.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
inline void MultiDArrayObject<DIM, VOXEL_TYPE, ARRAY_DIMENSION>::fillAll(VOXEL_TYPE voxelData) {
    vector<VOXEL_TYPE>::operator=(vector<VOXEL_TYPE>(this->size(), voxelData));
}

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
template<class INDEX_TYPE, class LOCAL_FUNCTION>
inline void MultiDArrayObject<DIM, VOXEL_TYPE, ARRAY_DIMENSION>::apply_voxelSlice(
    array<INDEX_TYPE, DIM> ijk, array<INDEX_TYPE, 2> limits, const LOCAL_FUNCTION& local_function) {
    static_assert(std::is_same_v<INDEX_TYPE, size_t> or std::is_same_v<INDEX_TYPE, long>);
    //////////////////
    ijk[DIM - 1] = 0;  // unimportant coordinate
    //////////////////
    if constexpr (std::is_same_v<INDEX_TYPE, long>) {  // beware of periodic replica
        array<size_t, DIM> ijkNew = vox::auxi::project_index_periodic<DIM>(ijk, this->arrayDim.getNbNodeBigGrid());
        vector<array<size_t, 2>> allSegments = vox::auxi::getSegmentsFromPeriodic(limits, this->arrayDim.getNbNodeBigGrid()[DIM - 1]);
        for (const auto& seg : allSegments) {
            apply_voxelSlice<size_t>(ijkNew, seg, local_function);
        }
    } else if constexpr (std::is_same_v<INDEX_TYPE, size_t>) {
        // in case of an extracted subgrid
        if constexpr (this->isSubArray) {
            if (not this->getArrayDimensions().doesCoverTorus()) {
                for (size_t i = 0; i < DIM - 1; i++) {
                    ijk[i] -= this->getArrayDimensions().getNMin()[i];
                }
                limits[0] -= this->getArrayDimensions().getNMin()[DIM - 1];
                limits[1] -= this->getArrayDimensions().getNMin()[DIM - 1];
            }
        }
        // in case of an extracted subgrid
        vox::auxi::loopOnIndexSlice<DIM>(ijk, limits, this->arrayDim.getNbNodeSubGrid(), [this, &local_function](size_t linearIndex) {
            this->apply_voxelLinear(linearIndex, local_function);
            });
    }
}

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
inline size_t MultiDArrayObject<DIM, VOXEL_TYPE, ARRAY_DIMENSION>::get_linear_index_periodic(const array<long, DIM>& ijk) const {
    return this->get_linear_index(vox::auxi::project_index_periodic<DIM>(ijk, this->arrayDim.getNbNodeBigGrid()));
}

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
inline size_t MultiDArrayObject<DIM, VOXEL_TYPE, ARRAY_DIMENSION>::get_linear_index(
    const array<size_t, DIM>& ijk) const {
    if constexpr (this->isSubArray) {
        if (this->getArrayDimensions().doesCoverTorus()) {
            return auxi::get_linear_index<DIM>(ijk, this->arrayDim.getNbNodeSubGrid());
        } else {
            auto ijkTilde = renormalize_index(ijk);
            return auxi::get_linear_index<DIM>(ijkTilde, this->arrayDim.getNbNodeSubGrid());
        }
    } else {
        return auxi::get_linear_index<DIM>(ijk, this->arrayDim.getNbNodeSubGrid());
    }
}

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
inline array<size_t, DIM> MultiDArrayObject<DIM, VOXEL_TYPE, ARRAY_DIMENSION>::renormalize_index(array<size_t, DIM> ijk) const {
    if constexpr (this->isSubArray) {
        for (size_t i = 0; i < DIM; i++) {
            ijk[i] -= this->getArrayDimensions().getNMin()[i];
        }
    }
    return ijk;
}

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
inline VOXEL_TYPE& MultiDArrayObject<DIM, VOXEL_TYPE, ARRAY_DIMENSION>::operator [](
    const array<size_t, DIM>& ijk) {
    return (*this)[this->get_linear_index(ijk)];
}

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
inline const VOXEL_TYPE& MultiDArrayObject<DIM, VOXEL_TYPE, ARRAY_DIMENSION>::operator [](
    const array<size_t, DIM>& ijk) const {
    return (*this)[this->get_linear_index(ijk)];
}

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
inline VOXEL_TYPE& MultiDArrayObject<DIM, VOXEL_TYPE, ARRAY_DIMENSION>::operator [](
    const array<long, DIM>& ijk) {
    return (*this)[this->get_linear_index_periodic(ijk)];
}

template<unsigned short DIM, class VOXEL_TYPE, class ARRAY_DIMENSION>
inline const VOXEL_TYPE& MultiDArrayObject<DIM, VOXEL_TYPE, ARRAY_DIMENSION>::operator [](
    const array<long, DIM>& ijk)  const {
    return (*this)[this->get_linear_index_periodic(ijk)];
}

// namespace auxi
template<unsigned short DIM, class FUNCTION>
inline void vox::auxi::loopOnIndexSlice(array<size_t, DIM> ijk, const array<size_t, 2>& limits, const array<size_t, DIM>& nbNodes,
    const FUNCTION& f) {
    static_assert(DIM == 2 or DIM == 3);
    size_t linearIndex = 0;
    //
    ijk[DIM - 1] = limits[0];
    linearIndex = auxi::get_linear_index<DIM>(ijk, nbNodes);
    //
    size_t limSlice_0 = limits[0];
    size_t limSlice_1 = limits[1];
    //
    for (size_t ijk_2 = limSlice_0; ijk_2 < limSlice_1; ijk_2++, linearIndex++) {
        f(linearIndex);
    }
}


template<unsigned short DIM, class FUNCTION>
inline void vox::auxi::loopOnIndex(const array<array<size_t, 2>, DIM>& limits, const array<size_t, DIM>& nbNodes,
    const FUNCTION& f) {
    static_assert(DIM == 2 or DIM == 3);
    if constexpr (DIM == 3) {
        size_t increment_0 = nbNodes[1] * nbNodes[2];
        size_t increment_1 = nbNodes[2];
        constexpr size_t increment_2 = 1;

        size_t ijk_0, ijk_1, ijk_2;
        size_t limits_0_0 = limits[0][0];
        size_t limits_0_1 = limits[0][1];
        size_t limits_1_0 = limits[1][0];
        size_t limits_1_1 = limits[1][1];
        size_t limits_2_0 = limits[2][0];
        size_t limits_2_1 = limits[2][1];

        array<size_t, DIM> ijk{ limits[0][0], limits[1][0], limits[2][0] };
        size_t linearIndex = auxi::get_linear_index<DIM>(ijk, nbNodes);
        for (ijk_0 = limits_0_0; ijk_0 < limits_0_1; ijk_0++, linearIndex += increment_0) {
            for (ijk_1 = limits_1_0; ijk_1 < limits_1_1; ijk_1++, linearIndex += increment_1) {
                for (ijk_2 = limits_2_0; ijk_2 < limits_2_1; ijk_2++, linearIndex += increment_2) {
                    f(linearIndex);
                }
                linearIndex -= (limits_2_1 - limits_2_0) * increment_2;
            }
            linearIndex -= (limits_1_1 - limits_1_0) * increment_1;
        }
    } else if constexpr (DIM == 2) {
        size_t increment_0 = nbNodes[1];
        constexpr size_t increment_1 = 1;
        size_t ijk_0, ijk_1;
        size_t limits_0_0 = limits[0][0];
        size_t limits_0_1 = limits[0][1];
        size_t limits_1_0 = limits[1][0];
        size_t limits_1_1 = limits[1][1];

        array<size_t, DIM> ijk{ limits[0][0], limits[1][0] };
        size_t linearIndex = auxi::get_linear_index<DIM>(ijk, nbNodes);
        for (ijk_0 = limits_0_0; ijk_0 < limits_0_1; ijk_0++, linearIndex += increment_0) {
            for (ijk_1 = limits_1_0; ijk_1 < limits_1_1; ijk_1++, linearIndex += increment_1) {
                f(linearIndex);
            }
            linearIndex -= (limits_1_1 - limits_1_0) * increment_1;
        }
    }
}



template<unsigned short DIM, class C>
inline array<size_t, DIM> vox::auxi::project_index_periodic(
    const array<long, DIM>& ijk, const C& nbNodes) {
    array<size_t, DIM> new_ijk{};
    for (size_t n = 0; n < DIM; n++) {
        new_ijk[n] = auxi_function::fast_modulo(ijk[n], nbNodes[n]);
    }
    return new_ijk;
}

template<unsigned short DIM, class C>
size_t auxi::get_linear_index_periodic(const array<long, DIM>& ijk, const C& nbNodes) {
    return get_linear_index<DIM>(project_index_periodic<DIM>(ijk, nbNodes), nbNodes);
}

template<unsigned short DIM, unsigned short INDEX, class C>
inline size_t get_linear_index_auxi_auxi(const array<size_t, DIM>& ijk,
    const C& nbNodes) {
    if constexpr (INDEX == 0) {
        return ijk[0];
    } else {
        return ijk[INDEX] + nbNodes[INDEX] * get_linear_index_auxi_auxi<DIM, INDEX - 1, C>(ijk, nbNodes);
    }
}

template<unsigned short DIM, class C>
size_t auxi::get_linear_index(const array<size_t, DIM>& ijk,
    const C& nbNodes) {
    assert(nbNodes.size() >= DIM);
    return get_linear_index_auxi_auxi<DIM, DIM - 1, C>(ijk, nbNodes);
}

template<unsigned short DIM>
inline vector<array<array<size_t, 2>, DIM> > auxi::getCuboidList(
    const array<array<long, 2>, DIM>& cuboidPerLimits,
    const array<size_t, DIM>& nbNodes) {
    static_assert(DIM == 2 or DIM == 3);
    ////////////////////////////////////
    array<vector<array<size_t, 2>>, DIM> allSegments{};
    array<size_t, DIM> numberOfSegments{};
    size_t numberOfCubes = 1;
    for (size_t i = 0; i < DIM; i++) {
        allSegments[i] = vox::auxi::getSegmentsFromPeriodic(cuboidPerLimits[i], nbNodes[i]);
        numberOfSegments[i] = allSegments[i].size();
        numberOfCubes *= numberOfSegments[i];
    }
    ////
    vector<array<array<size_t, 2>, DIM> > result(numberOfCubes);
    size_t counter = 0;
    ////
    loop<false, DIM>(numberOfSegments, [&result, &allSegments, &counter](const array<size_t, DIM> indexSegment) {
        for (size_t i = 0; i < DIM; i++) {
            result[counter][i] = allSegments[i][indexSegment[i]];
        }
        counter++;
        });
    return result;
}

}  // namespace vox
}  // namespace sac_de_billes

