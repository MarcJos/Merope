//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


namespace merope {
namespace homogenization {

template<Rule hr>
inline double homogenization::homog(const vector<double>& fracVol,
        const vector<double>& coeff) {
        if       constexpr (hr == Rule::Reuss)       return homogReuss(fracVol, coeff);
        else if  constexpr (hr == Rule::Voigt)       return homogVoigt(fracVol, coeff);
        else if constexpr (hr == Rule::Smallest)     return homogSmall(fracVol, coeff);
        else if constexpr (hr == Rule::Largest)      return homogLarge(fracVol, coeff);
        else  Merope_assert(false, "Invalid Homogenization Rule");
}

}  // namespace homogenization
}  // namespace merope



