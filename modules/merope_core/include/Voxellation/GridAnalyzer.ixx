//! Copyright : see license.txt
//!
//! \brief
#pragma once


namespace merope {
namespace vox {
namespace  grid_analyzer {

template<unsigned short DIM, class COMPOSITE, class>
std::map<PhaseType, double> compute_percentages(const CartesianGrid<DIM, COMPOSITE>& cartesianGrid) {
    std::map<PhaseType, double> res = {};
    double voxel_percentage = 1. / cartesianGrid.size();
    auto add_percentage = [&res](PhaseType phase, double volumeFraction) {
        if (res.find(phase) != res.end()) {
            res[phase] += volumeFraction;
        } else {
            res[phase] = volumeFraction;
        }
        };
    auto analyse_percentage = [add_percentage, voxel_percentage](auto composite_cell) {
        if constexpr (composite::is_Pure<COMPOSITE>) {
            add_percentage(composite_cell, 1. * voxel_percentage);
        } else if constexpr (composite::is_Iso<COMPOSITE> or composite::is_AnIso<COMPOSITE>) {
            for (const auto& fv : composite_cell) {
                add_percentage(fv.phase, fv.fracVol * voxel_percentage);
            }
        } else {
            Merope_assert(false, "Imposible");
        }
        };
    for (const auto& cell : cartesianGrid) {
        analyse_percentage(cell);
    }
    return res;
}

namespace auxi {
template<unsigned short DIM, class COMPOSITE>
std::map<PhaseType, double> compute_percentages_dynamical_error(const CartesianGrid<DIM, COMPOSITE>& cartesianGrid) {
    if constexpr (composite::is_composite<COMPOSITE>) {
        if constexpr (is_same_v<composite::Basic_Phase_Type<COMPOSITE>, PhaseType>) {
            return compute_percentages<DIM, COMPOSITE>(cartesianGrid);
        } else {
            throw runtime_error("Incorrect input type");
        }
    } else {
        throw runtime_error("Incorrect input type");
    }
}
}  // namespace aux

template<unsigned short DIM>
std::map<PhaseType, double> GridAnalyzer<DIM>::compute_percentages(const voxellizer::GridRepresentation<DIM>& gridRepresentation) const {
    std::map<PhaseType, double> res = {};
    auto compute_percentages_loc = [&res](const auto& grid) {
        res = auxi::compute_percentages_dynamical_error(grid);
        };
    gridRepresentation.try_on_all(compute_percentages_loc);
    return res;
}

string auxi::string_percentage(const std::map<PhaseType, double>& percentages) {
    string res = "";
    for (const auto& key_vec : percentages) {
        res += "phase " + std::to_string(get<0>(key_vec)) + " : " + std::to_string(100 * get<1>(key_vec)) + "%" + "\n";
    }
    return res;
}


}  // namespace  grid_analyzer
}  // namespace  vox
}  // namespace  merope