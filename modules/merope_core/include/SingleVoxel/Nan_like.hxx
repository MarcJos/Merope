//! Copyright : see license.txt
//!
//! \brief Functions for defining and recognizing nan = not a number phases.
//! Attempts to generalize for our context the notion of nan.

#pragma once

#include "../../../GenericMerope/StdHeaders.hxx"

namespace merope {

inline constexpr PhaseType magic_nan_long(){
    return numeric_limits<PhaseType>::min();
}

template<class BasicType>
BasicType make_nan_like() {
    if constexpr (std::is_same_v<BasicType, PhaseType>) {
        return magic_nan_long();
    } else if constexpr (std::is_same_v<BasicType, double>) {
        return std::nan("");
    } else {
        Merope_static_error(BasicType, "Impossible");
    }
}

template<class BasicType>
bool is_nan_like(BasicType number) {
    if constexpr (std::is_same_v<BasicType, PhaseType>) {
        return number == magic_nan_long();
    } else if constexpr (std::is_same_v<BasicType, double>) {
        return std::isnan(number);
    } else {
        Merope_static_error(BasicType, "Impossible");
    }
}

}  // namespace  merope
