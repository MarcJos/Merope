//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"
#include "../Voxellation/Voxellizer.hxx"

namespace merope {
namespace vox {

//! @warning : not programmed yet

namespace grid_analyzer {

template<unsigned short DIM, class COMPOSITE,
    class = std::enable_if_t<is_same_v<composite::Basic_Phase_Type<COMPOSITE>, PhaseType>>>
std::map<PhaseType, double>  compute_percentages(const CartesianGrid<DIM, COMPOSITE>& cartesianGrid);

namespace auxi {
inline string string_percentage(const std::map<PhaseType, double>& percentages);
}  // namespace  auxi

template<unsigned short DIM>
class GridAnalyzer {
public:
    void print_percentages(const voxellizer::GridRepresentation<DIM>& gridRepresentation) const {
        cout << auxi::string_percentage(compute_percentages(gridRepresentation));
    }

    std::map<PhaseType, double> compute_percentages(const voxellizer::GridRepresentation<DIM>& gridRepresentation) const;
};

}  // namespace  grid_analyzer

}  // namespace  vox
}  // namespace  merope

#include "../Voxellation/GridAnalyzer.ixx"