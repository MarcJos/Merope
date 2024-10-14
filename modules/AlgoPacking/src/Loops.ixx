//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace sac_de_billes {
using namespace std;

/*!
 * \file   exemple.cxx
 * \brief
 * \author th202608
 * \date   25/03/2021
 */

 /*template<array<long, 3> I_BEGIN, array<long, 3> I_END>
  constexpr size_t getSizeIndices(){
  size_t sizeIndices = 1;
  for(short i=0; i<3; i++){
  sizeIndices *= I_END[i] - I_BEGIN[i];
  }
  return sizeIndices;
  }


  template<array<long, 3>* I_BEGIN, array<long, 3>* I_END>
  constexpr array<array<size_t, 3u>, getSizeIndices<I_BEGIN, I_END>()> getIndices() {
  auto c = size_t { };
  auto indices = std::array<std::array<size_t, 3u>, getSizeIndices<I_BEGIN, I_END>()> { };
  for (size_t i1 = (*I_BEGIN)[0]; i1 != (*I_END)[0]; ++i1) {
  for (size_t i2 = (*I_BEGIN)[1]; i2 != (*I_END)[1]; ++i2) {
  for (size_t i3 = (*I_BEGIN)[2]; i3 != (*I_END)[2]; ++i3, ++c) {
  indices[c] = { i1, i2, i3 };
  }
  }
  }
  return indices;
  }*/

template<unsigned short I, unsigned short DIM, class INDEX_TYPE, class C>
void putNewIndices(const C& limits, size_t& counter, array<INDEX_TYPE, DIM>& index, vector<array<INDEX_TYPE, DIM>>& result) {
    for (index[I] = limits[I][0]; index[I] < limits[I][1]; index[I]++) {
        if constexpr (I == 0) {
            result[counter] = index;
            counter++;
        } else {
            putNewIndices<I - 1, DIM, INDEX_TYPE, C>(limits, counter, index, result);
        }
    }
}

template<unsigned short DIM1, unsigned short DIM2, class INDEX_TYPE, class C>
vector<array<INDEX_TYPE, DIM1>> getAllIndices_multi(C limits) {
    assert(limits.size() >= DIM2);
    size_t resultSize = 1;
    for (size_t i = 0; i < DIM2; i++) {
        resultSize *= (limits[i][1] - limits[i][0]);
    }
    vector<array<INDEX_TYPE, DIM1>> result{ };
    result.resize(resultSize);
    size_t counter = 0;
    array<INDEX_TYPE, DIM1> index;
    putNewIndices<DIM2 - 1, DIM1, INDEX_TYPE, C>(limits, counter, index, result);
    return result;
}

template<unsigned short DIM, class INDEX_TYPE, class C>
vector<array<INDEX_TYPE, DIM>> getAllIndices(C limits) {
    return getAllIndices_multi<DIM, DIM, INDEX_TYPE, C>(limits);
}

template<unsigned short DIM1, unsigned short DIM2, class INDEX_TYPE, class C>
vector<array<INDEX_TYPE, DIM1>> getAllIndices_from0_multi(C limits) {
    using arrElemType = decltype(+limits[0]);
    array<array<arrElemType, 2>, limits.size()> newLimits{};
    for (size_t i = 0; i < limits.size(); i++) {
        newLimits[i][0] = 0;
        newLimits[i][1] = limits[i];
    }
    return getAllIndices_multi<DIM1, DIM2, INDEX_TYPE>(newLimits);
}

template<unsigned short DIM, class INDEX_TYPE, class C>
vector<array<INDEX_TYPE, DIM>> getAllIndices_from0(C limits) {
    return getAllIndices_from0_multi<DIM, DIM, INDEX_TYPE>(limits);
}

template<size_t I1, size_t I2, size_t I3>
constexpr std::array<std::array<size_t, 3u>, I1* I2* I3> getIndices() {
    auto c = size_t{ };
    auto indices = std::array<std::array<size_t, 3u>, I1* I2* I3> { };
    for (size_t i1 = 0; i1 != I1; ++i1) {
        for (size_t i2 = 0; i2 != I2; ++i2) {
            for (size_t i3 = 0; i3 != I3; ++i3, ++c) {
                indices[c] = { i1, i2, i3 };
            }
        }
    }
    return indices;
}

template<size_t I1, size_t I2>
constexpr std::array<std::array<size_t, 2u>, I1* I2> getIndices() {
    auto c = size_t{ };
    auto indices = std::array<std::array<size_t, 2u>, I1* I2> { };
    for (size_t i1 = 0; i1 != I1; ++i1) {
        for (size_t i2 = 0; i2 != I2; ++i2, ++c) {
            indices[c] = { i1, i2 };
        }
    }
    return indices;
}

template<unsigned short DIM, class C, class Ctilde>
inline vector<array<C, DIM>> cartesianProduct(Ctilde indicesLimits) {
    // implicitly have Ctilde = array<vector<C>, DIM2>
    assert(indicesLimits.size() >= DIM);
    ////////////////////////////
    vector<array<C, DIM>> result{};
    if constexpr (DIM == 1) {
        for (const auto& elem : indicesLimits[0]) {
            result.push_back(array<C, DIM>{elem});
        }
    } else {
        vector<array<C, DIM - 1>> previousProduct = cartesianProduct<DIM - 1, C>(indicesLimits);
        array<C, DIM> currentElem;
        for (const array<C, DIM - 1>&elem2 : previousProduct) {
            for (const C& elem : indicesLimits[DIM - 1]) {
                //
                for (size_t i = 0; i < DIM - 1;i++) {
                    currentElem[i] = elem2[i];
                }
                currentElem[DIM - 1] = elem;
                //
                result.push_back(currentElem);
            }
        }
    }
    return result;
}
}  // namespace sac_de_billes


