//! Copyright : see license.txt
//!
//! \brief
//
#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Grid/GridManipulations.hxx"


#include "MeropeNamespace.hxx"


namespace merope {
namespace vox {

array<vector<double>, 2> gridAuxi::getTabCoeff(
    const VoxelPhaseFrac& phaseFrac,
    const vector<double>& coeff) {
    size_t nbPhase = phaseFrac.size();
    vector<double> tab_coeffs(nbPhase, 0);
    vector<double> tab_fracVols(nbPhase, 0);
    for (size_t i = 0; i < nbPhase; i++) {
        tab_coeffs[i] = coeff[phaseFrac[i].phase];
        tab_fracVols[i] = phaseFrac[i].fracVol;
    }
    return array<vector<double>, 2> { tab_fracVols, tab_coeffs };
}

array<double, 2> gridAuxi::translateMask(
    const VoxelPhaseFrac& mask_i) {
    array<double, 2> result{ 0., 0. };
    for (const auto& phfv : mask_i) {
        if (phfv.phase == 1) {
            result[1] += phfv.fracVol;
        }
        else {
            result[0] += phfv.fracVol;
        }
    }
    return result;
}

} //namespace vox
} // namespace merope



