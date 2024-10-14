//! Copyright : see license.txt
//!
//! \brief
//

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Field/CovSum.hxx"


#include "MeropeNamespace.hxx"


namespace merope {
namespace gaussianField {

// Default constructor of a Covariance Function
// r: Range
// s: Sill
CovIso::CovIso(const double r_i, const double s_i) :
    r(r_i), s(s_i) {}

// Default constructor of a Covariance Function
CovIso::~CovIso() {}

// 1D Covariance Function
// hx: Distance
double CovIso::cov(const double hx) const {
    return s * covIso(fabs(hx / r));
}

// 2D Covariance Function
// hx,hy: Distance
double CovIso::cov(const double hx, const double hy) const {
    return s * covIso(hypot(hx, hy) / r);
}

// 3D Covariance Function
// hx,hy,hz: Distance
double CovIso::cov(const double hx, const double hy, const double hz) const {
    return s * covIso(sqrt(hx * hx + hy * hy + hz * hz) / r);
}

// Default constructor of the exponential covariance
// r: Range
// s: Sill
CovExpo::CovExpo(const double r_i, const double s_i) :
    CovIso(r_i, s_i) {}

// Isotrope covariance function
// h: Distance
inline double CovExpo::covIso(const double h) const {
    if (h > 700) return 0;
    return exp(-h);
}

// Default constructor of the spherical covariance
// r: Range
// s: Sill
CovSpherical::CovSpherical(const double r_i, const double s_i) :
    CovIso(r_i, s_i) {}

// Isotrope covariance function
// h: Distance
inline double CovSpherical::covIso(const double h) const {
    if (h > 1) return 0;
    else return 1 - (1.5 - 0.5 * h * h) * h;
}

// Default constructor of the spherical covariance
// r: Range
// s: Sill
CovGaussian::CovGaussian(const double r_i, const double s_i) :
    CovIso(r_i, s_i) {}

// Isotrope covariance function
// h: Distance
inline double CovGaussian::covIso(const double h) const {
    if (h > 26) return 0;
    return exp(-h * h);
}

void CovSum::add(CovType covType, vector<double> modelParameters) {
    switch (covType) {
    case CovType::Exponential:
    {
        this->add(new CovExpo(modelParameters[0], modelParameters[1]));
        break;
    }
    case CovType::Gaussian:
    {
        this->add(new CovGaussian(modelParameters[0], modelParameters[1]));
        break;
    }
    case CovType::Spherical:
    {
        this->add(new CovSpherical(modelParameters[0], modelParameters[1]));
        break;
    }
    default:
        throw(invalid_argument("CovSum_add: model unknown"));
    }
}

// Add a convariance function
// c: New convariance function
void CovSum::add(CovIso* c) {
    push_back(c);
}

// 1D Covariance Function
// hx: Distance
double CovSum::cov(const double hx) const {
    double sum = 0;
    for (CovSum::const_iterator i = begin(); i != end(); ++i) {
        sum += (*i)->cov(hx);
    }
    return sum;
}

// 2D Covariance Function
// hx,hy: Distance
double CovSum::cov(const double hx, const double hy) const {
    double sum = 0;
    for (CovSum::const_iterator i = begin(); i != end(); ++i) {
        sum += (*i)->cov(hx, hy);
    }
    return sum;
}

// 3D Covariance Function
// hx,hy,hz: Distance
double CovSum::cov(const double hx, const double hy, const double hz) const {
    double sum = 0;
    for (CovSum::const_iterator i = begin(); i != end(); ++i) {
        sum += (*i)->cov(hx, hy, hz);
    }
    return sum;
}

// Default Destructor of a covariance sum
CovSum::~CovSum() {
    for (CovSum::iterator i = begin(); i != end(); ++i) {
        delete* i;
    }
}


}  // namespace gaussianField
}  // namespace merope

