//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../MeropeNamespace.hxx"


namespace merope {
namespace homogenization {

template<Rule hr>
inline double homogenization::homog(const vector<double>& fracVol,
        const vector<double>& coeff) {
        if       constexpr (hr == Rule::Reuss)       return homogReuss(fracVol, coeff);
        else if  constexpr (hr == Rule::Voigt)       return homogVoigt(fracVol, coeff);
        else if constexpr (hr == Rule::Smallest)     return homogSmall(fracVol, coeff);
        else if constexpr (hr == Rule::Largest)      return homogLarge(fracVol, coeff);
        else  throw invalid_argument(__PRETTY_FUNCTION__);
}

}  // namespace homogenization
}  // namespace merope



