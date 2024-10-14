//! Copyright : see license.txt
//!
//! \briefIncludes the other loops
//
#pragma once

#include "StdHeaders.hxx"

namespace sac_de_billes {
using namespace std;

//! \return all the indices ([index[0], index[1], ..., index[DIM]]) such that limits[j][0] <= i[j] < limits[j][1]
template<unsigned short DIM1, unsigned short DIM2, class INDEX_TYPE, class C>
vector<array<INDEX_TYPE, DIM1>> getAllIndices_multi(C limits);
//! case DIM1 = DIM2
template<unsigned short DIM, class INDEX_TYPE, class C>
vector<array<INDEX_TYPE, DIM>> getAllIndices(C limits);

//! \return all the indices ([index[0], index[1], ..., index[DIM]]) such that 0 <= i[j] < limits[j][1]
template<unsigned short DIM1, unsigned short DIM2, class INDEX_TYPE, class C>
vector<array<INDEX_TYPE, DIM1>> getAllIndices_from0_multi(C limits);
//! case DIM1 = DIM2
template<unsigned short DIM, class INDEX_TYPE, class C>
vector<array<INDEX_TYPE, DIM>> getAllIndices_from0(C limits);

//! return a cartesian product of segments
// implicitly have Ctilde = array<vector<C>, DIM2>
template<unsigned short DIM, class C, class Ctilde>
vector<array<C, DIM>> cartesianProduct(Ctilde indicesLimits);

}  // namespace sac_de_billes

#include "Loops.ixx"
#include "Loops_2.ixx"


