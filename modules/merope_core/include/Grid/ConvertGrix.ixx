//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_CONVERTGRIX_IXX_
#define GRID_CONVERTGRIX_IXX_


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM, class C2, class FUNCTION>
inline void vox::convertGrid::apply_inplace(
    CartesianGrid<DIM, C2>& grid0, FUNCTION rule) {
#pragma omp parallel for firstprivate(rule)
    for (size_t i = 0; i < grid0.size(); i++) {
        rule(grid0[i]);
    }
}

template<unsigned short DIM, class C1, class C2, class FUNCTION>
inline CartesianGrid<DIM, C1> vox::convertGrid::auxi::localConvertCartesian(const CartesianGrid<DIM, C2>& grid0, FUNCTION rule) {
    CartesianGrid<DIM, C1> resultGrid(grid0.getGridParameters());
#pragma omp parallel for firstprivate(rule)
    for (size_t i = 0; i < grid0.size(); i++) {
        resultGrid[i] = rule(grid0[i]);
    }
    return resultGrid;
}

template<unsigned short DIM, class C1, class C2, class FUNCTION>
inline CartesianGrid<DIM, C1> vox::convertGrid::localConvert(const CartesianGrid<DIM, C2>& grid0, const FUNCTION& rule) {
    return auxi::localConvertCartesian<DIM, C1, C2, FUNCTION>(grid0, rule);
}

template<unsigned short DIM, class PHASE_TYPE>
CartesianGrid<DIM, gridAuxi::ListPhaseFrac<PHASE_TYPE>> vox::convertGrid::fromPhaseToFracVol(const CartesianGrid<DIM, PHASE_TYPE>& grid0) {
    return localConvert<DIM, gridAuxi::ListPhaseFrac<PHASE_TYPE>, PHASE_TYPE>(grid0, [](const PHASE_TYPE& phase) {
        return (gridAuxi::ListPhaseFrac<PHASE_TYPE>(
            { auxi_SphereCollection::PhaseFrac<PHASE_TYPE>(phase, 1.) }));
        });
}

template<unsigned short DIM>
inline vector<vector<tuple<vox::VTK_PHASE, double>>> vox::convertGrid::fromCartesianToVector(
    CartesianGrid<DIM, VoxelPhaseFrac> grid0) {
    return static_cast<vector<vector<tuple<vox::VTK_PHASE, double>>>>(localConvert<DIM, vector<tuple<vox::VTK_PHASE, double>>, VoxelPhaseFrac>(grid0,
        [](const auto& voxelPhaseFrac) {
            vector<tuple<vox::VTK_PHASE, double>> result = {};
            for (const auto& phaseFrac : voxelPhaseFrac) {
                result.push_back(make_tuple(phaseFrac.phase, phaseFrac.fracVol));
            }
            return result;
        }));
}

template<unsigned short DIM>
inline CartesianGrid<DIM, VTK_PHASE> vox::convertGrid::fromFieldToPhase(const CartesianGrid<DIM, double>& gridField, vector<double>& fieldValues) {
    double coeffmax = *max_element(gridField.begin(), gridField.end());
    double coeffmin = *min_element(gridField.begin(), gridField.end());
    double stepCoeff = (coeffmax - coeffmin) / NB_PHASE_USHORT;
    // Put all phases in buckets
    CartesianGrid<DIM, VTK_PHASE> result = localConvert<DIM, VTK_PHASE, double>(gridField,
        [coeffmin, coeffmax, stepCoeff](const auto& fieldValue) {
            VTK_PHASE phase = rint((fieldValue - coeffmin) / stepCoeff);
            if (phase > vox::convertGrid::NB_PHASE_USHORT) {
                phase = vox::convertGrid::NB_PHASE_USHORT;
            }
            return phase;
        });
    // Modify accordingly the values of the coefficient
    fieldValues = vector<double>(NB_PHASE_USHORT + 1);
    for (size_t i = 0; i < fieldValues.size(); i++) {
        fieldValues[i] = coeffmin + i * stepCoeff;
    }
    // manipulate the result such that there only remains the non-empty buckets
    convertGrid::renormalizeWithCoefficients(result, fieldValues);
    //
    return result;
}

template<unsigned short DIM>
inline void vox::convertGrid::renormalizeWithCoefficients(CartesianGrid<DIM, VTK_PHASE>& grid0, vector<double>& fieldValues) {
    unique_ptr<TabPhaseCoeff<VTK_PHASE>> res;
    res.reset(new TabPhaseCoeff<VTK_PHASE>(grid0, fieldValues));
    res->verifyCoherent();
    res->removeUnusedPhase();
}

} //namespace vox
} // namespace merope

#endif /* GRID_CONVERTGRIX_IXX_ */
