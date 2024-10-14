//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"
#include "../MeropeNamespace.hxx"
#include "../Grid/CompositeVoxelTypes.hxx"


namespace merope {
namespace vox {
namespace gridAuxi {

template<class PHASE_FRAC_>
class ListPhaseFrac : public vector<PHASE_FRAC_> {
    //! Class used to represent on a given point the different volume fractions of each phase
    //! Ideally, these volume fractions sum up to 100%
public:
    using PHASE_FRAC = PHASE_FRAC_;
    static_assert(is_PhaseFrac<PHASE_FRAC> or is_PhaseFracNormal<PHASE_FRAC>);
    //! constructor
    ListPhaseFrac() = default;
    //! constructor
    ListPhaseFrac(PHASE_FRAC::PHASE_TYPE phase) : vector<PHASE_FRAC>(
        { PHASE_FRAC(phase, 1.) }) {}
    //! constructor
    ListPhaseFrac(initializer_list<PHASE_FRAC> lst)
        : vector<PHASE_FRAC>(lst) {}
    //! adds another phasefrac, but maintaining minimal memory requirements
    void add(const PHASE_FRAC& pf) { this->push_back(pf); }
    //! merge phases that are present twice (adding the volume fractions)
    void merge(double merge_criterion = 1e-6);
    //! renormalizes the volume fractions to 100%, eventually add the matrixPhase
    //! \param matrixPhase : phase to be eventually added to get 1 volume fraction
    //! \param is_there_matrix : if true, there is a matrix
    void renormalize(bool is_there_matrix = true, PHASE_FRAC::PHASE_TYPE matrixPhase = 0);
private:
    //! @brief renormalizes the volume fractions to 100%
    //! by dividing all volume fractions by the given inverseTotalVolumeFraction
    void renormalize_by_multiply(double inverseTotalVolumeFraction);
};

}  // namespace gridAuxi

namespace composite {

//! @return a global vector representative of the anisotropy of the composite voxel
//! @param aniso_composite : list of (phase, volume_fraction, normal)
template<class PHASE_FRAC>
Point<PHASE_FRAC::DIM> compute_global_normal(const gridAuxi::ListPhaseFrac<PHASE_FRAC>& aniso_composite);
}  // namespace  composite

}  // namespace vox
}  // namespace merope

#include "../Grid/ListPhaseFrac.ixx"


