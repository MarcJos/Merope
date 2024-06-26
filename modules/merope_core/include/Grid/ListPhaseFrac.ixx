//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_LISTPHASEFRAC_IXX_
#define GRID_LISTPHASEFRAC_IXX_

#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {
namespace gridAuxi {

template<class TYPE_PHASE>
inline void ListPhaseFrac<TYPE_PHASE>::merge(double merge_criterion) {
    if (this->size() < 2) return; // size = 0, 1
    ///////////////////////////
    sort(this->begin(), this->end());
    size_t i_current = 0;
    for (size_t i_next = 1; i_next < this->size(); i_next++) {
        if (std::abs((*this)[i_current].phase - (*this)[i_next].phase) < merge_criterion) {
            (*this)[i_current].fracVol += (*this)[i_next].fracVol;
        } else {
            i_current++;
            (*this)[i_current] = (*this)[i_next];
        }
    }
    this->resize(i_current + 1);
}

template<class TYPE_PHASE>
inline void  ListPhaseFrac<TYPE_PHASE>::renormalize(bool is_there_matrix, TYPE_PHASE matrixPhase) {
    if ((*this).size() == 0) {
        if (not(is_there_matrix)) {
            throw runtime_error("Unexpected : no matrix phase has been set, and the voxel is empty");
        } else {
            (*this).push_back(
                auxi_SphereCollection::PhaseFrac<TYPE_PHASE> { matrixPhase, 1. });
        }
        return;
    }
    double totalVolumeFraction = 0.;
    for (const auto& phfv : (*this)) {
        totalVolumeFraction += phfv.fracVol;
    }
    if (abs(totalVolumeFraction - 1) < geomTools::EPSILON) { // totalVolumeFraction close to 1
        return;
    }
    double inverseTotalVolumeFraction = 1. / totalVolumeFraction;
    if (not(is_there_matrix) or totalVolumeFraction > 1) {
        renormalize_by_multiply(inverseTotalVolumeFraction);
    } else { // totalVolumeFraction too small
        (*this).push_back(auxi_SphereCollection::PhaseFrac<TYPE_PHASE> { matrixPhase,
            1. - totalVolumeFraction});
    }
}

template<class TYPE_PHASE>
inline void  ListPhaseFrac<TYPE_PHASE>::renormalize_by_multiply(double inverseTotalVolumeFraction) {
    for (auto& phfv : (*this)) {
        phfv.fracVol *= inverseTotalVolumeFraction;
    }
}

} // namespace vox::gridAuxi
} // namespace vox
} // namespace merope


#endif /* GRID_LISTPHASEFRAC_IXX_ */
