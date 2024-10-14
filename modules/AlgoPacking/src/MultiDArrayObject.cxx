//! Copyright : see license.txt
//!
//! \brief
#include "StdHeaders.hxx"
#include "MultiDArrayObject.hxx"

namespace merope {

vector<array<size_t, 2>> vox::auxi::getSegmentsFromPeriodic(const array<long, 2>& perLimits, size_t nbNode) {
    assert(perLimits[0] <= perLimits[1]);
    if (perLimits[1] - perLimits[0] >= nbNode) {
        return vector<array<size_t, 2>>{array<size_t, 2>{0, nbNode}};
    } else {
        size_t k = 0;  // for perLimits[0] - k * nbNode
        long i0_long = perLimits[0];
        while (i0_long < 0) {
            k--;
            i0_long += nbNode;
        }
        while (i0_long >= nbNode) {
            k++;
            i0_long -= nbNode;
        }
        size_t i0 = static_cast<size_t>(i0_long);
        size_t i1 = static_cast<size_t>(perLimits[1] - k * nbNode);
        if (i1 <= nbNode) {
            return vector<array<size_t, 2>>{array<size_t, 2>{i0, i1}};
        } else {
            return vector<array<size_t, 2>>{array<size_t, 2>{i0, nbNode}, array<size_t, 2>{0, i1 - nbNode}};
        }
    }
}

}  // namespace merope
