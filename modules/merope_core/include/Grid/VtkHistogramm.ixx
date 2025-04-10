//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


namespace merope {
namespace vox {

template<class PHASE_TYPE, class COEFF_TYPE>
bool TabPhaseCoeff<PHASE_TYPE, COEFF_TYPE>::verifyCoherent() const {
    PHASE_TYPE maxPhase = *max_element(this->phases.begin(), this->phases.end());
    bool isCoherent = maxPhase < this->coeff.size();
    if (not isCoherent) {
        cerr << "Maximal phase measured " << maxPhase << endl;
        cerr << "Size of the coefficient table" << this->coeff.size() << endl;
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("TabPhaseCoeff not coherent!");
    }
    return isCoherent;
}

template<class PHASE_TYPE, class COEFF_TYPE>
void TabPhaseCoeff<PHASE_TYPE, COEFF_TYPE>::guaranteeIncreasingValues() {
    bool needSorting = not std::is_sorted(coeff.begin(), coeff.end());
    if (not needSorting) {
        return;
    }
    vector<size_t> indexCoeff(coeff.size());
    for (size_t i = 0; i < indexCoeff.size(); i++) {
        indexCoeff[i] = i;
    }
    sort(indexCoeff.begin(), indexCoeff.end(),
        [this](const auto& i, const auto& j) {
            return this->coeff[i] < this->coeff[j];
        }
    );
    // update coeff
    auto newCoeff = coeff;
    for (size_t i = 0; i < coeff.size(); i++) {
        newCoeff[indexCoeff[i]] = coeff[i];
    }
    coeff = newCoeff;
    // update phases
    for (auto& ph : phases) {
        ph = indexCoeff[ph];
    }
}

template<class PHASE_TYPE, class COEFF_TYPE>
void TabPhaseCoeff<PHASE_TYPE, COEFF_TYPE>::removeUnusedPhase() {
    assert(verifyCoherent());
    guaranteeIncreasingValues();
    ///////////////////////////////////
    vector<size_t> phaseList(coeff.size());
    vector<size_t> voxelPerPhases(coeff.size(), 0);
    for (size_t i = 0; i < phaseList.size(); i++) {
        phaseList[i] = i;
    }
    // sees which phase are present
    for (auto& phase : phases) {
        (voxelPerPhases[phase])++;
    }
    // Remove the empy phases
    size_t number_of_void_phases = 0;
    for (size_t i = 0; i < phaseList.size(); i++) {
        if (voxelPerPhases[i] == 0) {
            number_of_void_phases++;
        }
        phaseList[i] -= number_of_void_phases;
    }  // at the end, thePhase[phase] links phase with \tilde{phase}, which such that any \tilde{phase} is at least represented once in tabPhaseCoeff
    size_t number_of_phases = phaseList[phaseList.size() - 1] + 1;
    // set the new phase
    for (auto& onePhase : phases) {
        onePhase = phaseList[onePhase];
    }
    // set the new coefficient
    vector<COEFF_TYPE> newCoeff(number_of_phases);
    size_t index = 0;
    for (size_t i = 0; i < phaseList.size(); i++) {
        if (voxelPerPhases[i] > 0) {
            newCoeff[index] = coeff[i];
            index++;
        }
    }
    coeff = newCoeff;
}

}  // namespace vox
}  // namespace merope



