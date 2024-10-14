//! Copyright : see license.txt
//!
//! \brief
//!


#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Geometry/Hyper.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

Hyper::Hyper(const double xa, const double ya, const double xb, const double yb,
    const double Ra, const double Rb, const double EPS2_i) :
    EPS2(EPS2_i) {
    // Point milieu
    xc = (xa + xb) / 2;
    yc = (ya + yb) / 2;
    double dx = xb - xa;
    double dy = yb - ya;
    double d2AB = dx * dx + dy * dy;
    double dAB = sqrt(d2AB);
    if (Ra > Rb) {
        DRs2 = (Ra - Rb) / 2;
        ux = dx / dAB;
        uy = dy / dAB;
    } else {
        DRs2 = (Rb - Ra) / 2;
        ux = -dx / dAB;
        uy = -dy / dAB;
    }
    m2 = d2AB / (4 * DRs2 * DRs2) - 1;
    m = sqrt(m2);
    m4 = m2 * m2;
}

Hyper::Hyper(const double xa, const double ya, const double za, const double xb,
    const double yb, const double zb, const double Ra, const double Rb,
    const double EPS2_i) :
    EPS2(EPS2_i) {
    // Point milieu
    xc = (xa + xb) / 2;
    yc = (ya + yb) / 2;
    zc = (za + zb) / 2;
    double dx = xb - xa;
    double dy = yb - ya;
    double dz = zb - za;
    double d2AB = dx * dx + dy * dy + dz * dz;
    double dAB = sqrt(d2AB);
    if (Ra > Rb) {
        DRs2 = (Ra - Rb) / 2;
        ux = dx / dAB;
        uy = dy / dAB;
        uz = dz / dAB;
    } else {
        DRs2 = (Rb - Ra) / 2;
        ux = -dx / dAB;
        uy = -dy / dAB;
        uz = -dz / dAB;
    }
    m2 = d2AB / (4 * DRs2 * DRs2) - 1;
    m = sqrt(m2);
    m4 = m2 * m2;
}

void Hyper::projM(const double um, const double vm, double& up,
    double& vp) const {
    double d0 = m * vm + um;
    vp = m * (m2 * d0 - sqrt(d0 * d0 + m4 - 1)) / (m4 - 1);
    up = um + m * (vm - vp);
}

void Hyper::projD(const double um, const double vm, double& up,
    double& vp) const {
    double d0 = -m * vm + um;
    vp = -m * (m2 * d0 - sqrt(d0 * d0 + m4 - 1)) / (m4 - 1);
    up = um - m * (vm - vp);
}

void Hyper::init(const double um, const double vm, double& up,
    double& vp) const {
    // Domaine 1
    if ((um - 1 + m * vm > 0) and (vm > 0)) {
        projM(um, vm, up, vp);
        return;
    }
    // Domaine 2
    if ((um - 1 - m * vm > 0) and (vm <= 0)) {
        projD(um, vm, up, vp);
        return;
    }
    // Domaine 3
    up = 1;
    vp = 0;
}

double Hyper::oneIter(const double um, const double vm, double& u,
    double& v) const {
    double P = m2 * (1 - u * u) + v * v;
    double Q = (u - um) * v + m2 * (v - vm) * u;
    // double Q=2*((u-um)*v+m2*(v-vm)*u);
    double Pu = -2 * m2 * u;
    double Pv = 2 * v;
    double Qu = v + m2 * (v - vm);
    double Qv = u - um + m2 * u;
    // double Qu=2*(v+m2*(v-vm));
    // double Qv=2*(u-um+m2*u);
    double D = Pu * Qv - Pv * Qu;

    double du = (Qv * P - Pv * Q) / D;
    u -= du;
    double dv = (Qu * P - Pu * Q) / D;
    v += dv;
    return du * du + dv * dv;
}

double Hyper::distance(const double um, const double vm) const {
    double up, vp;
    // Initialisation
    init(um, vm, up, vp);
    // Iterations de Newton
    unsigned char i = 0;
    double d2;
    do {
        d2 = oneIter(um, vm, up, vp);
        if (++i > 10)
            throw(runtime_error("Hyper::distance: pas de convergence"));
    } while (d2 > EPS2);
    up -= um;
    vp -= vm;
    return sqrt(up * up + vp * vp);
}

double Hyper::distanceRO(double xm, double ym) const {
    // Recentrage et normalisation
    xm = (xm - xc) / DRs2;
    ym = (ym - yc) / DRs2;
    // Coordonnees du point (xm,ym) dans le repere u,v où v est perpendiculaire a u calcule dans le constructeur -> v(-uy,ux)
    double um = ux * xm + uy * ym;
    double vm = -uy * xm + ux * ym;
    return DRs2 * distance(um, vm);
}

double Hyper::distanceRO(double xm, double ym, double zm) const {
    // Recentrage et normalisation
    xm = (xm - xc) / DRs2;
    ym = (ym - yc) / DRs2;
    zm = (zm - zc) / DRs2;

    // Longueur CM
    double d2CM = xm * xm + ym * ym + zm * zm;
    // double wx = uy*zm - ym*uz;
    // double wy = uz*xm - zm*ux;
    // double wz = ux*ym - xm*uy;
    // double dw = sqrt(wx*wx + wy*wy + wz*wz);
    // wx /= dw;
    // wy /= dw;
    // wz /= dw;

    // // Vecteur v = w vectoriel u
    // double vx = wy*uz - uy*wz;
    // double vy = wz*ux - uz*wx;
    // double vz = wx*uy - ux*wy;

    // Coordonnees du point (xm,ym,zm) dans le repere u,v,w où wm = 0, um = CM.u, vm = CM.v ou vm se calcul a partir de d2CM = um*um + vm*vm
    double um = ux * xm + uy * ym + uz * zm;
    double vm = sqrt(d2CM - um * um);
    // double vm_ = vx*xm + vy*ym + vz*zm;
    // double wm = wx*xm + wy*ym + wz*zm;

    return DRs2 * distance(um, vm);
}
}  // namespace merope
