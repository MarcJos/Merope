//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../Grid/GridTypesBase.hxx"


namespace merope {
namespace vox {
namespace composite {

template<unsigned short DIM, class PHASE_OUT, class COMPOSITE,
    class CONVERSION, typename = enable_if_t<is_composite<COMPOSITE>>>
Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> apply_conversion_local(const COMPOSITE& local_data,
    const CONVERSION& conversion) {
    if constexpr (is_Pure<COMPOSITE>) {
        return conversion(local_data);
    }
    //
    Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> result{};
    if constexpr (is_Iso<COMPOSITE>) {
        for (const auto& pf : local_data) {
            result.push_back(PhaseFrac<PHASE_OUT>(conversion(pf.phase), pf.fracVol));
        }
    }
    if constexpr (is_AnIso<COMPOSITE>) {
        for (const auto& pf : local_data) {
            result.push_back(PhaseFracNormal<DIM, PHASE_OUT>(conversion(pf.phase), pf.fracVol, pf.normal));
        }
    }
    return result;
}

template<unsigned short DIM, class COMPOSITE, class TEXTURER,
    class PHASE_OUT = double, typename = enable_if_t<is_composite<COMPOSITE>>>
Change_Type_Composite<DIM, COMPOSITE, PHASE_OUT> apply_texture_loc(const COMPOSITE& local_data,
    const TEXTURER& texturer, const Point<DIM>& x) {
    auto conversion = [&](const auto& phase) {
        return texturer(x, phase);
        };
    return apply_conversion_local<DIM, PHASE_OUT>(local_data, conversion);
}

}  // namespace  composite
}  // namespace vox
}  // namespace merope

#include "GridTypesBase.ixx"


