//! Copyright : see license.txt
//!
//! \brief Geostatistic data needed to simulate a periodic medium
//!


#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Grid/Geostat.hxx"
#include "../../AlgoPacking/src/Geometry/GeomTypes.hxx"

#include "MeropeNamespace.hxx"


namespace merope {
namespace gaussianField {

template<typename M>
void StepDis::FromModel(const M& mod, const double x0, const double dx,
    const double eps) {
    double F1 = 0;
    double xp = x0;
    xi.push_back(xp);
    while (F1 < (1 - eps)) {
        // New step
        xp += dx;
        xi.push_back(xp);
        // Cumulative distribution function
        F1 = mod.CDF(xp);
        Fi.push_back(F1);
    }
    xp += dx;
    xi.push_back(xp);
    F1 = 1;
    Fi.push_back(F1);
}


// Probability density function of a Normal Random Variable
// x: Value
inline double f_g(const double x) {
    static const double us2p = 1 / sqrt(2 * m_PI);
    double a = -x * x * 0.5;
    if (a < -300) return 0;
    return exp(a) * us2p;
}

// Cumulative distribution function of a Normal Random Variable
// x: Value
inline double F_g(const double g) {
    constexpr double my_M_SQRT1_2 = 1 / sqrt(2);
    return 0.5 * (erf(g * my_M_SQRT1_2) + 1);
}

// RV

RV::~RV() {}

double RV::inverseA(const double x) const {
    return inverseF(F_g(x));
}

// StepDis

StepDis::StepDis() {}

void StepDis::read(string nom) {
    ifstream fin(nom);
    string line;
    double DF0 = -1, F = 0, x0 = 0;
    do {
        getline(fin, line);
        if ('#' != line[0]) {
            double x, DF;
            if (2 == sscanf(&line[0], "%lf%lf", &x, &DF)) {
                // Previous record
                if (DF0 >= 0) {
                    xi.push_back(x0);
                    F += DF0;
                    Fi.push_back(F);
                }
                DF0 = DF;
                x0 = x;
            }
        }
    } while (!fin.eof());
    xi.push_back(x0);
    // Last value of the CDF must be equal to 1
    Fi.back() = 1;
}

StepDis::StepDis(string nom) {
    read(nom);
}

size_t StepDis::size() const {
    return xi.size();
}

void StepDis::clear() {
    xi.clear();
    Fi.clear();
}

double StepDis::inverseF(double Fv) const {
    vector<double>::const_iterator x = xi.begin(), F = Fi.begin();
    if (Fv < 0) return *x;
    if (Fv > 1) Fv = 1;

    double xO = *(x++);
    double FO = 0;
    double xN = *(x++);
    double FN = *(F++);

    while (F != Fi.end()) {
        if (FN > Fv) break;
        xO = xN;
        FO = FN;
        xN = *(x++);
        FN = *(F++);
    }
    return xO + (xN - xO) * (Fv - FO) / (FN - FO);
}

void StepDis::setApprox() const {}

// GaussRV

GaussRV::GaussRV(const double m_i, const double s_i, const double w_i):
    m(m_i), s(s_i), w(w_i) {}

GaussRV::~GaussRV() {}

double GaussRV::PDF(const double x) const {
    return w * f_g((x - m) / s) / s;
}

double GaussRV::CDF(const double x) const {
    return w * F_g((x - m) / s);
}

inline void TGauss::MAJ() {
    Fl = F_g((l - m) / s);
    A = 1 / (1 - Fl);
    Fl *= w;
}

TGauss::TGauss(const double m_i, const double s_i, const double w_i):
    GaussRV(m_i, s_i, w_i), l(0) {
    MAJ();
}

double TGauss::PDF(const double x) const {
    if (x < l) {
        return 0;
    }
    return A * GaussRV::PDF(x);
}

double TGauss::CDF(const double x) const {
    if (x < l) {
        return 0.;
    }
    else {
        return A * (GaussRV::CDF(x) - Fl);
    }
}

LogNorm::LogNorm(const double m_i, const double s_i, const double w_i):
    GaussRV(m_i, s_i, w_i) {}

double LogNorm::PDF(const double x) const {
    if (x <= 0) {
        return 0;
    }
    return GaussRV::PDF(log(x)) / x;

}

double LogNorm::CDF(const double x) const {
    if (x <= 0) {
        return 0.;
    }
    else {
        return GaussRV::CDF(log(x));
    }
}

SofERV::~SofERV() {
    // Désallouer la mémoire
    for (auto& i : *this) {
        delete i;
    }
    vector<GaussRV*>::clear();
}

void SofERV::add(GaussRV* e) {
    push_back(e);
    sd.clear();
}

double SofERV::PDF(const double x) const {
    double s = 0;
    for (const auto& i : *this) {
        s += i->PDF(x);
    }
    return s;
}

double SofERV::CDF(const double x) const {
    double s = 0;
    for (const auto& i : *this) {
        s += i->CDF(x);
    }
    return s;
}

void SofERV::setApprox() const {
    if (!sd.size()) {
        sd.FromModel(*this, 0, 0.02);
    }
}

double SofERV::inverseF(const double u) const {
    if (u <= 0) return 0;
    if (u >= 1) throw(logic_error("SofERV::inverseF: u must be < 1"));

    double t1 = sd.inverseF(u);
    for (unsigned short i = 0; i < 20; ++i) {
        double t0 = t1;
        t1 = t0 + (u - CDF(t0)) / PDF(t0);
        if (fabs(t1 - t0) < 1e-6) {
            return t1;
        }
    }
    throw(runtime_error("SofERV::inverseF: No convergence"));
    return 0;
}

} // namespace gaussianField
} // namespace merope

