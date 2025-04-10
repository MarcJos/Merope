//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../Geometry/include/BasicGeometricOperations.hxx"
#include "../../../GenericTools/CPP_Functions.hxx"


namespace merope {
namespace vox {
namespace gridAuxi {
using namespace sac_de_billes;

template<class PHASE_FRAC>
void ListPhaseFrac<PHASE_FRAC>::postProcess(
    const MatrixPhaseHolder<typename PHASE_FRAC::PHASE_TYPE>& matrixPhaseHolder) {
    this->merge();
    this->renormalize(matrixPhaseHolder.is_there_matrix(), matrixPhaseHolder.getMatrixPhase());
}


template<class PHASE_FRAC>
void ListPhaseFrac<PHASE_FRAC>::merge(double merge_criterion) {
    if (this->size() < 2) return;  // size = 0, 1
    ///////////////////////////
    auto merge2phases = [](auto& currentPhase, const auto& otherPhase) {
        if constexpr (is_PhaseFracNormal<PHASE_FRAC>) {
            RenormPoint<PHASE_FRAC::DIM> new_normal(
                currentPhase.fracVol * currentPhase.normal
                + otherPhase.fracVol * otherPhase.normal);
            currentPhase.normal = new_normal.getPoint();
        }
        currentPhase.fracVol += otherPhase.fracVol;
        };

    sort(this->begin(), this->end());
    size_t i_current = 0;
    for (size_t i_next = 1; i_next < this->size(); i_next++) {
        if (abs((*this)[i_current].phase - (*this)[i_next].phase) < merge_criterion) {
            merge2phases((*this)[i_current], (*this)[i_next]);
        } else {
            i_current++;
            (*this)[i_current] = (*this)[i_next];
        }
    }
    this->resize(i_current + 1);
}

template<class PHASE_FRAC>
void  ListPhaseFrac<PHASE_FRAC>::renormalize(bool is_there_matrix, typename ListPhaseFrac<PHASE_FRAC>::PHASE_FRAC::PHASE_TYPE matrixPhase) {
    // case empty voxel
    if (this->size() == 0) {
        Merope_assert(is_there_matrix,
            "Unexpected : no matrix phase has been set, and the voxel is empty");
        this->push_back(PHASE_FRAC{ matrixPhase, 1. });
        return;
    }
    // case not-full voxel
    double totalVolumeFraction = 0.;
    for (const auto& phfv : (*this)) {
        totalVolumeFraction += phfv.fracVol;
    }
    if (abs(totalVolumeFraction - 1) < geomTools::EPSILON) {  // totalVolumeFraction close to 1
        return;
    }
    double inverseTotalVolumeFraction = 1. / totalVolumeFraction;
    // should fill in the voxel
    if (not(is_there_matrix) or totalVolumeFraction > 1) {
        renormalize_by_multiply(inverseTotalVolumeFraction);
    } else {  // totalVolumeFraction too small
        if constexpr (is_PhaseFrac<PHASE_FRAC>) {
            this->push_back(PHASE_FRAC{ matrixPhase, 1. - totalVolumeFraction });
        } else if constexpr (is_PhaseFracNormal<PHASE_FRAC>) {
            auto average_vector = vox::composite::compute_global_normal(*this);
            this->push_back(PHASE_FRAC(matrixPhase, 1. - totalVolumeFraction, average_vector));
        }
    }
}

template<class PHASE_TYPE>
void  ListPhaseFrac<PHASE_TYPE>::renormalize_by_multiply(double inverseTotalVolumeFraction) {
    for (auto& phfv : (*this)) {
        phfv.fracVol *= inverseTotalVolumeFraction;
    }
}

}  // namespace vox::gridAuxi

template<class PHASE_FRAC>
Point<PHASE_FRAC::DIM> composite::compute_global_normal(const gridAuxi::ListPhaseFrac<PHASE_FRAC>& aniso_composite) {
    constexpr unsigned short DIM = PHASE_FRAC::DIM;
    Point<DIM> average_vector = create_array<DIM>(0.);
    for (const auto& current_pfn : aniso_composite) {
        int sign_normal = (geomTools::prodScal<DIM>(current_pfn.normal, average_vector) > 0) ? 1 : -1;
        average_vector += sign_normal * current_pfn.fracVol * current_pfn.normal;
    }
    geomTools::renormalize<DIM>(average_vector);
    return average_vector;
}

}  // namespace vox
}  // namespace merope



