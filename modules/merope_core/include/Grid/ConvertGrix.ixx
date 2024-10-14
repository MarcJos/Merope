//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


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

template<bool use_open_mp, unsigned short DIM, class C1, class C2, class FUNCTION>
inline CartesianGrid<DIM, C1> vox::convertGrid::auxi::localConvertCartesian(const CartesianGrid<DIM, C2>& grid0, FUNCTION rule) {
    CartesianGrid<DIM, C1> resultGrid(grid0.getGridParameters());
    if constexpr (use_open_mp) {
#pragma omp parallel for firstprivate(rule)
        for (size_t i = 0; i < grid0.size(); i++) {
            resultGrid[i] = rule(grid0[i]);
        }
    } else {
        for (size_t i = 0; i < grid0.size(); i++) {
            resultGrid[i] = rule(grid0[i]);
        }
    }
    return resultGrid;
}

template<bool use_open_mp, unsigned short DIM, class C1, class C2, class FUNCTION>
inline CartesianGrid<DIM, C1> vox::convertGrid::localConvert(const CartesianGrid<DIM, C2>& grid0, const FUNCTION& rule) {
    return auxi::localConvertCartesian<use_open_mp, DIM, C1, C2, FUNCTION>(grid0, rule);
}

template<unsigned short DIM, class C1, class C2, class TEXTURER>
CartesianGrid<DIM, C1> convertGrid::applyFunctionDependingOnX(const CartesianGrid<DIM, C2>& grid0, const TEXTURER& texturer) {
    vox::CartesianGrid<DIM, C1> result(grid0.getGridParameters());
    sac_de_billes::loop<true>(grid0.getArrayDimensions().getNbNodes(),
        [&grid0, &result, &texturer](const auto& ijk) {
            result[ijk] = texturer(grid0.getCenterVoxel(ijk),
                grid0[ijk]);
        });
    return result;
}

template<unsigned short DIM, class PHASE_TYPE>
CartesianGrid<DIM, composite::Iso<PHASE_TYPE>> vox::convertGrid::fromPureToIso(const CartesianGrid<DIM, PHASE_TYPE>& grid0) {
    return localConvert<true, DIM, composite::Iso<PHASE_TYPE>, PHASE_TYPE>(grid0, [](const PHASE_TYPE& phase) {
        return composite::Iso<PHASE_TYPE>(phase);
        });
}

template<unsigned short DIM, class PHASE_TYPE>
CartesianGrid<DIM, vox::composite::AnIso<DIM, PHASE_TYPE>> vox::convertGrid::fromPureToAnIso(const CartesianGrid<DIM, PHASE_TYPE>& grid0) {
    return localConvert<DIM, gridAuxi::ListPhaseFrac<PHASE_TYPE>, PHASE_TYPE>(grid0, [](const PHASE_TYPE& phase) {
        return vox::composite::AnIso<DIM, PHASE_TYPE>(phase);
        });
}

template<unsigned short DIM, class COMPOSITE>
CartesianGrid<DIM, composite::to_stl_format<COMPOSITE>> vox::convertGrid::convert_to_stl_format(const CartesianGrid<DIM, COMPOSITE>& grid0) {
    return localConvert<true, DIM, composite::to_stl_format<COMPOSITE>, COMPOSITE>(grid0,
        &composite::convert_to_stl_format<COMPOSITE, DIM>);
}

template<unsigned short DIM, class VOXEL_TYPE>
vector<VOXEL_TYPE> vox::convertGrid::linearize(const CartesianGrid<DIM, VOXEL_TYPE>& grid0) {
    return static_cast<vector<VOXEL_TYPE>>(grid0);
}

template<unsigned short DIM>
CartesianGrid<DIM, PhaseType> vox::convertGrid::fromFieldToPhase(const CartesianGrid<DIM, double>& gridField, vector<double>& fieldValues) {
    double coeffmax = *max_element(gridField.begin(), gridField.end()) + 1000 * numeric_limits<double>::min();
    double coeffmin = *min_element(gridField.begin(), gridField.end()) - 1000 * numeric_limits<double>::min();
    double stepCoeff = (coeffmax - coeffmin) / NB_PHASE_USHORT + 1000 * numeric_limits<double>::min();
    // Put all phases in buckets
    CartesianGrid<DIM, PhaseType> result = localConvert<true, DIM, PhaseType, double>(gridField,
        [coeffmin, coeffmax, stepCoeff](const auto& fieldValue) {
            PhaseType phase = rint((fieldValue - coeffmin) / stepCoeff);
            // checks
            {
                if (phase > vox::convertGrid::NB_PHASE_USHORT) {
                    phase = vox::convertGrid::NB_PHASE_USHORT;
                }
                if (phase < 0) {
                    cerr << __PRETTY_FUNCTION__ << endl;
                    cerr << "fieldValue : " << fieldValue << " vs " << "coeffmin : " << coeffmin << " , " << "stepCoeff : " << stepCoeff << endl;
                }
            }
            //
            return phase;
        });
    // Modify accordingly the values of the coefficient
    PhaseType max_phase = *(std::max_element(result.begin(), result.end()));
    fieldValues = vector<double>(max_phase + 1);
    for (size_t i = 0; i < fieldValues.size(); i++) {
        fieldValues[i] = coeffmin + i * stepCoeff;
    }
    // manipulate the result such that there only remains the non-empty buckets
    convertGrid::renormalizeWithCoefficients(result, fieldValues);
    //
    return result;
}

template<unsigned short DIM, class COEFF_TYPE, class C>
void vox::convertGrid::renormalizeWithCoefficients(CartesianGrid<DIM, PhaseType>& grid0, vector<COEFF_TYPE>& fieldValues) {
    TabPhaseCoeff<PhaseType, COEFF_TYPE> res(grid0, fieldValues);
    res.verifyCoherent();
    res.removeUnusedPhase();
}

template<unsigned short DIM>
void  vox::convertGrid::removeUnusedPhase(CartesianGrid<DIM, PhaseType>& grid, vector<PhaseType>& phaseValues) {
    auto max_elem = *max_element(grid.begin(), grid.end());
    phaseValues.resize(max_elem + 1);
    for (size_t i = 0; i < phaseValues.size(); i++) {
        phaseValues[i] = i;
    }
    convertGrid::renormalizeWithCoefficients(grid, phaseValues);
}

}  // namespace vox
}  // namespace merope


