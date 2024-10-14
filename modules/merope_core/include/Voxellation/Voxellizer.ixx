//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../Grid/ManipOn_GridTypeBases.hxx"

namespace merope {
namespace vox {
namespace voxellizer {

template<homogenization::Rule HOMOG_RULE, unsigned short DIM>
GridField<DIM> applyHomogRule_T(const CartesianGrid<DIM, composite::Iso<double>>& phaseFracVol) {
    CartesianGrid<DIM, double> gridField = convertGrid::localConvert<true, DIM, double, composite::Iso<double>>
        (phaseFracVol,
            [](const auto& valueFracVolLoc) {
                array<vector<double>, 2> input;
                for (size_t i = 0; i < valueFracVolLoc.size(); i++) {
                    input[0].push_back(valueFracVolLoc[i].fracVol);
                    input[1].push_back(valueFracVolLoc[i].phase);
                }
                return homogenization::homog<HOMOG_RULE>(input[0], input[1]);
            }
        );
    return gridField;
}

template<homogenization::Rule HOMOG_RULE, unsigned short DIM>
GridField<DIM> applyHomogRule_T(const CartesianGrid<DIM, composite::Iso<PhaseType>>& phaseFracVol,
    const vector<double>& pureCoeffs) {
    CartesianGrid<DIM, double> gridField = convertGrid::localConvert<true, DIM, double, composite::Iso<PhaseType>>(phaseFracVol,
        [&pureCoeffs](const auto& phaseFracVolLoc) {
            auto input = gridAuxi::getTabCoeff(phaseFracVolLoc, pureCoeffs);
            return homogenization::homog<HOMOG_RULE>(input[0], input[1]);
        }
    );
    return gridField;
}

template<unsigned short DIM, class COMPOSITE, class TEXTURER, typename T>
CartesianGrid<DIM, composite::Change_Type_Composite<DIM, COMPOSITE, double>> apply_texture(const CartesianGrid<DIM, COMPOSITE>& grid0, const TEXTURER& texturer) {
    return convertGrid::applyFunctionDependingOnX<DIM, composite::Change_Type_Composite<DIM, COMPOSITE, double>>(grid0,
        [&texturer](const auto& x, const auto& local_data) {
            return composite::apply_texture_loc<DIM, COMPOSITE>(local_data, texturer, x);
        });
}


template<unsigned short DIM, class COMPOSITE>
CartesianGrid<DIM, composite::Change_Type_Composite<DIM, COMPOSITE, double>> apply_coefficients(const CartesianGrid<DIM, COMPOSITE>& grid, const vector<double>& pureCoeffs) {
    if constexpr (composite::is_Pure<COMPOSITE>) {
        auto conversion = [&pureCoeffs](auto i) {
            return pureCoeffs[i];
            };
        return convertGrid::auxi::localConvertCartesian<true, DIM, double>(grid, conversion);
    } else if constexpr (composite::is_Iso<COMPOSITE>) {
        auto conversion = [&pureCoeffs](const auto& listphaseFrac) {
            composite::Iso<double> res{};
            res.reserve(listphaseFrac.size());
            for (const auto& pfn : listphaseFrac) {
                res.push_back(PhaseFrac<double>(pureCoeffs[pfn.phase], pfn.fracVol));
            }
            return res;
            };
        return convertGrid::auxi::localConvertCartesian<true, DIM, composite::Iso<double>>(grid, conversion);
    } else if constexpr (composite::is_AnIso<COMPOSITE>) {
        auto conversion = [&pureCoeffs](const auto& listphaseFrac) {
            composite::AnIso<DIM, double> res{};
            res.reserve(listphaseFrac.size());
            for (const auto& pfn : listphaseFrac) {
                res.push_back(PhaseFracNormal<DIM, double>(pureCoeffs[pfn.phase], pfn.fracVol, pfn.normal));
            }
            return res;
            };
        return convertGrid::auxi::localConvertCartesian<true, DIM, composite::AnIso<DIM, double>>(grid, conversion);
    } else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected");
    }
}


}  // namespace  voxellizer
} /* namespace vox */
}  // namespace merope