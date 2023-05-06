//! Copyright : see license.txt
//!
//! \brief

#ifndef SRC_ELLIPSE_CXX
#define SRC_ELLIPSE_CXX

#include "Geometry/Ellipse.hxx"
#include "Geometry/Area.hxx"
#include "AmbiantSpace.hxx"

double sac_de_billes::ellipseAux::computeChordArea(const Ellipse<2>& ellipse, Point<2> x_0, Point<2> x_1) {
    RenormPoint<2> n_x(ellipse.axes[0]);
    RenormPoint<2> n_y(ellipse.axes[1]);
    x_0 -= ellipse.center;
    x_1 -= ellipse.center;
    double theta0 = std::atan2(geomTools::prodScal<2>(n_y, x_0), geomTools::prodScal<2>(n_x, x_0));
    double theta1 = std::atan2(geomTools::prodScal<2>(n_y, x_1), geomTools::prodScal<2>(n_x, x_0));
    // Renormalization in case of incoherent angles
    if (theta1 < theta0) {
        theta1 += 2 * M_PI;
    }
    double result = computeAngularArea(ellipse, theta1) - computeAngularArea(ellipse, theta0); // only angular part
    // decide whether triangle should be added or removed
    if (theta1 - theta0 < M_PI) {
        result -= geomTools::area::triangle<2>({ 0, 0 }, x_0, x_1);
    } else {
        result += geomTools::area::triangle<2>({ 0, 0 }, x_0, x_1);
    }
    //
    return result;
}

double sac_de_billes::ellipseAux::computeAngularArea(const Ellipse<2>& ellipse, double theta) {
    //! see https://www.geometrictools.com/Documentation/AreaIntersectingEllipses.pdf
    double a = ellipse.getAlphas()[0];
    double b = ellipse.getAlphas()[1];
    return 0.5 * a * b * (theta -
        atan((b - a) * sin(2 * theta)
            / ((b + a) + (b - a) * cos(2 * theta))));
}

#endif // SRC_ELLIPSE_CXX