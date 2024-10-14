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
    const composite::Iso<PhaseType>& phaseFrac,
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
    const composite::Iso<PhaseType>& mask_i) {
    array<double, 2> result{ 0., 0. };
    for (const auto& phfv : mask_i) {
        if (phfv.phase == 1) {
            result[1] += phfv.fracVol;
        } else {
            result[0] += phfv.fracVol;
        }
    }
    return result;
}

array<double, 2> gridAuxi::translateMask(const composite::Iso<double>& voxelMask) {
    // test preconditions
    double density_error = 1e-6;
    if ((voxelMask.size() != 1) or (abs(voxelMask[0].fracVol - 1) > density_error)
        or (voxelMask[0].phase < -density_error) or (voxelMask[0].phase - 1 > density_error)) {
        cerr << std::boolalpha << endl;
        cerr << (voxelMask.size() != 1) << endl;
        cerr << (abs(voxelMask[0].fracVol - 1) > density_error) << endl;
        cerr << (voxelMask[0].phase < 0) << endl;
        cerr << (voxelMask[0].phase > 1) << endl;
        cerr << "Voxel Mask phase " << voxelMask[0].phase << endl;
        cerr << "Voxel Mask phase - 1 " << voxelMask[0].phase - 1 << endl;
        cerr << "Voxel Mask fracVol " << voxelMask[0].fracVol << endl;
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected. Mask should represent a density");
    }

    double phi = max(min(1., voxelMask[0].phase), 0.);
    array<double, 2> result{ 1 - phi, phi };
    return result;
}

}  // namespace vox
}  // namespace merope



