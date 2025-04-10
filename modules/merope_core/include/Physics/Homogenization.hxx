//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"
#include "../../../GenericTools/CPP_Functions.hxx"


namespace merope {
//! namespace containing information relative to thermal homogenization rules
namespace homogenization {
enum class Rule {
    Reuss,      // harmonic mean
    Voigt,      // arithmetic mean
    Smallest,   // take the smallest coefficient in the cell
    Largest     // take the largest coefficient in the cell
};

//! see Rule
double homogReuss(const vector<double>& fracVol, const vector<double>& coeff);
//! see Rule
double homogVoigt(const vector<double>& fracVol, const vector<double>& coeff);
//! see Rule
double homogSmall(const vector<double>& fracVol, const vector<double>& coeff);
//! see Rule
double homogLarge(const vector<double>& fracVol, const vector<double>& coeff);
//! chooses between the rules. Compile-time.
template<Rule hr>
double homog(const vector<double>& fracVol, const vector<double>& coeff);
//! chooses between the rules. Not compile-time.
std::function<double(const vector<double>&, const vector<double>&)> getHomogRule(Rule rule);

}  // namespace homogenization
}  // namespace merope


#include "../Physics/Homogenization.ixx"


