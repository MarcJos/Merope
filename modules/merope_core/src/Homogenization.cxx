//! Copyright : see license.txt
//!
//! \brief
//

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Physics/Homogenization.hxx"

#include "MeropeNamespace.hxx"


namespace merope {

namespace homogenization {

double homogReuss(const vector<double>& fracVol,
    const vector<double>& coeff) {
    double coeffTot = 0;
    for (size_t i = 0; i < fracVol.size(); i++) {
        coeffTot += fracVol[i] * 1. / coeff[i];
    }
    return 1. / coeffTot;
}

double homogVoigt(const vector<double>& fracVol,
    const vector<double>& coeff) {
    double coeffTot = 0;
    for (size_t i = 0; i < fracVol.size(); i++) {
        coeffTot += fracVol[i] * coeff[i];
    }
    return coeffTot;
}

double homogSmall(const vector<double>&,
    const vector<double>& coeff) {
    return *(min_element(coeff.begin(), coeff.end()));
}

double homogLarge(const vector<double>&,
    const vector<double>& coeff) {
    return *(max_element(coeff.begin(), coeff.end()));
}

std::function<double(const vector<double>&, const vector<double>&)> getHomogRule(homogenization::Rule hr) {
    if (hr == Rule::Reuss)        return [](const vector<double>& fracVol, const vector<double>& coeff) {return homogReuss(fracVol, coeff);};
    else if (hr == Rule::Voigt)        return [](const vector<double>& fracVol, const vector<double>& coeff) {return homogVoigt(fracVol, coeff);};
    else if (hr == Rule::Smallest)     return [](const vector<double>& fracVol, const vector<double>& coeff) {return homogSmall(fracVol, coeff);};
    else if (hr == Rule::Largest)      return [](const vector<double>& fracVol, const vector<double>& coeff) {return homogLarge(fracVol, coeff);};
    else  throw invalid_argument(__PRETTY_FUNCTION__);
}

}  // namespace homogenization
}  // namespace merope
