//! Copyright : see license.txt
//!
//! \brief

#pragma once

namespace merope {
namespace vox {
namespace composite {

template<class COMPOSITE, unsigned short DIM>
to_stl_format<COMPOSITE> convert_to_stl_format(const COMPOSITE& composite_) {
    static_assert(is_Pure<COMPOSITE> or is_Iso<COMPOSITE> or is_AnIso<COMPOSITE>);
    using PHASE_TYPE = Basic_Phase_Type<COMPOSITE>;
    if constexpr (is_Pure<COMPOSITE>) {
        return composite_;
    }
    if constexpr (is_Iso<COMPOSITE>) {
        composite::stl_format_Iso<PHASE_TYPE> result = {};
        for (const auto& phaseFrac : composite_) {
            result.push_back(make_tuple(phaseFrac.phase, phaseFrac.fracVol));
        }
        return result;
    }
    if constexpr (is_AnIso<COMPOSITE>) {
        composite::stl_format_Iso<PHASE_TYPE> result_phaseFrac = {};
        for (const auto& phaseFrac : composite_) {
            result_phaseFrac.push_back(make_tuple(phaseFrac.phase, phaseFrac.fracVol));
        }
        Point<DIM> normal = vox::composite::compute_global_normal(composite_);
        composite::stl_format_AnIso<DIM, PHASE_TYPE> result = { result_phaseFrac, normal };
        return result;
    }
}


}  // namespace  composite
}  // namespace vox
}  // namespace merope