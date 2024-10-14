//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<class PHASE_TYPE, class COEFF_TYPE>
class TabPhaseCoeff {
    //! auxiliary class to make sure that the final grid has phases [0, N-1], with at least one voxel per phase
public:
    TabPhaseCoeff(vector<PHASE_TYPE>& phases_, vector<COEFF_TYPE>& coeff_) :
        phases{ phases_ }, coeff{ coeff_ } {}
    //! verifies if coeff_ is sufficiently large
    bool verifyCoherent() const;
    //! at the end, each coeff is present at least in one phase[i]
    void removeUnusedPhase();
private:
    //! sort coeff in increasing order, modify accordingly the phases
    void guaranteeIncreasingValues();
    //! represents the phase grid
    vector<PHASE_TYPE>& phases;
    //! the values of the coeff at voxel is coeff[phases[i]]
    vector<COEFF_TYPE>& coeff;
};

}  // namespace vox
}  // namespace merope

#include "../Grid/VtkHistogramm.ixx"


