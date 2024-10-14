//! Copyright : see license.txt
//!
//! \brief
//!
//!


#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Geometry/Rotation.hxx"
#include "../../AlgoPacking/src/Geometry/GeomTypes.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

inline void Rotation2D::VecProd() {
    a[2] = -a[1];  // a[1,2]
    a[3] = a[0];  // a[2,2]
}

Rotation2D::Rotation2D() {
    a[0] = a[1] = a[2] = a[3] = NAN;
}

Rotation2D::Rotation2D(const CastemReal* const v, const bool verif) {
    if (verif) {
        // Normalisation
        CastemReal n = sqrt(v[0] * v[0] + v[1] * v[1]);
        a[0] = v[0] / n;  // a[1,1]
        a[1] = v[1] / n;  // a[2,1]
    } else {
        a[0] = v[0];  // a[1,1]
        a[1] = v[1];  // a[2,1]
    }
    // Second vector: vectorial product
    VecProd();
}

Rotation2D::Rotation2D(const double theta) {
    EulerToAxis(theta);
}

void Rotation2D::randomEuler(default_random_engine& eng) {
    uniform_real_distribution<> dist(0, m_PI);
    double theta = dist(eng);

    EulerToAxis(theta);
}

void Rotation2D::EulerToAxis(const double theta) {
    // First vector
    a[0] = cos(theta);  // a[1,1]
    a[1] = sin(theta);  // a[2,1]
    VecProd();  // Second vector: vectorial product
}

const CastemReal* Rotation2D::getMat() const {
    return a;
}

void Rotation2D::direct(const CastemReal* const vm,
    CastemReal* const vg) const {
    vg[0] = a[0] * vm[0] + a[2] * vm[1];
    vg[1] = a[1] * vm[0] + a[3] * vm[1];
}

void Rotation2D::reverse(const CastemReal* const vg,
    CastemReal* const vm) const {
    vm[0] = a[0] * vg[0] + a[1] * vg[1];
    vm[1] = a[2] * vg[0] + a[3] * vg[1];
}

void Rotation2D::print(ofstream& f) const {
    f << setprecision(auxiRotations::PRINTED_DOUBLE_PRECISION);
    f << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << endl;
}

Rotation3D::Rotation3D() {
    a[0] = a[1] = a[2] = a[3] = a[4] = a[5] = a[6] = a[7] = a[8] = NAN;
}

inline void Rotation3D::VecProd() {
    // Third vector
    a[6] = a[1] * a[5] - a[4] * a[2];  // a[1,3]
    a[7] = a[2] * a[3] - a[5] * a[0];  // a[2,3]
    a[8] = a[0] * a[4] - a[3] * a[1];  // a[3,3]
}

Rotation3D::Rotation3D(const CastemReal* const v, const bool verif) {
    if (verif) {
        // First vector
        // Normalisation
        CastemReal n = sqrt(v[0] * v[0] + v[1] * v[1] + v[2] * v[2]);
        a[0] = v[0] / n;  // a[1,1]
        a[1] = v[1] / n;  // a[2,1]
        a[2] = v[2] / n;  // a[3,1]

        // Third vector
        a[6] = v[1] * v[5] - v[2] * v[4];  // a[1,3]
        a[7] = v[2] * v[3] - v[0] * v[5];  // a[2,3]
        a[8] = v[0] * v[4] - v[1] * v[3];  // a[3,3]
        // Normalisation
        n = sqrt(a[6] * a[6] + a[7] * a[7] + a[8] * a[8]);
        a[6] /= n;
        a[7] /= n;
        a[8] /= n;

        // Second vector
        a[3] = a[7] * a[2] - a[8] * a[1];  // a[1,2]
        a[4] = a[8] * a[0] - a[6] * a[2];  // a[2,2]
        a[5] = a[6] * a[1] - a[7] * a[0];  // a[3,2]
    } else {
        a[0] = v[0];  // a[1,1]
        a[1] = v[1];  // a[2,1]
        a[2] = v[2];  // a[3,1]

        // Second vector
        a[3] = v[3];  // a[1,2]
        a[4] = v[4];  // a[2,2]
        a[5] = v[5];  // a[3,2]

        // Third vector: vectorial product
        VecProd();
    }
}

const CastemReal* Rotation3D::getMat() const {
    return a;
}

Rotation3D::Rotation3D(const double phi, const double theta, const double psi) {
    EulerToAxis(phi, theta, psi);
}

void Rotation3D::randomEuler(default_random_engine& eng) {
    uniform_real_distribution<> dist(0, 2 * m_PI), dist2(0, 1);
    // Équivalence Bunge-Euler : phi_1 = psi, phi = theta, phi_2 = phi
    double psi = dist(eng);
    double theta = acos(1 - 2 * dist2(eng));
    double phi = dist(eng);

    EulerToAxis(phi, theta, psi);
}

void Rotation3D::EulerToAxis(const double phi, const double theta,
    const double psi) {
    const double sPhi = sin(phi);
    const double cPhi = cos(phi);
    const double sTheta = sin(theta);
    const double cTheta = cos(theta);
    const double sPsi = sin(psi);
    const double cPsi = cos(psi);

    // First vector
    a[0] = cPsi * cPhi - sPsi * cTheta * sPhi;  // a[1,1]
    a[1] = cPsi * sPhi + sPsi * cTheta * cPhi;  // a[2,1]
    a[2] = sPsi * sTheta;  // a[3,1]

    // Second vector
    a[3] = -sPsi * cPhi - cPsi * cTheta * sPhi;  // a[1,2]
    a[4] = -sPsi * sPhi + cPsi * cTheta * cPhi;  // a[2,2]
    a[5] = cPsi * sTheta;  // a[3,2]

    // Third vector: vectorial product
    VecProd();
}

void Rotation3D::direct(const CastemReal* const vm,
    CastemReal* const vg) const {
    vg[0] = a[0] * vm[0] + a[3] * vm[1] + a[6] * vm[2];
    vg[1] = a[1] * vm[0] + a[4] * vm[1] + a[7] * vm[2];
    vg[2] = a[2] * vm[0] + a[5] * vm[1] + a[8] * vm[2];
}

void Rotation3D::reverse(const CastemReal* const vg,
    CastemReal* const vm) const {
    vm[0] = a[0] * vg[0] + a[1] * vg[1] + a[2] * vg[2];
    vm[1] = a[3] * vg[0] + a[4] * vg[1] + a[5] * vg[2];
    vm[2] = a[6] * vg[0] + a[7] * vg[1] + a[8] * vg[2];
}

void Rotation3D::print(ofstream& f) const {
    f << setprecision(auxiRotations::PRINTED_DOUBLE_PRECISION);
    f << a[0] << " " << a[1] << " " << a[2] << " " << a[3] << " " << a[4] << " "
        << a[5] << endl;
}

}  // namespace merope

