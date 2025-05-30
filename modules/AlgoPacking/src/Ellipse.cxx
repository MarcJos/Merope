//! Copyright : see license.txt
//!
//! \brief

#include "../../Geometry/include/Ellipse.hxx"
#include "../../Geometry/include/Area.hxx"
#include "../../Geometry/include/AmbiantSpace.hxx"

namespace merope {

double ellipseAux::angleOnEllipse(const Ellipse<2>& ellipse, Point<2> pointOnEllipse) {
    RenormPoint<2> n_x(ellipse.axes[0]);
    RenormPoint<2> n_y(ellipse.axes[1]);
    pointOnEllipse -= ellipse.center;
    return std::atan2(geomTools::prodScal<2>(n_y, pointOnEllipse), geomTools::prodScal<2>(n_x, pointOnEllipse));
}

double ellipseAux::computeChordArea(const Ellipse<2>& ellipse, Point<2> x_0, Point<2> x_1) {
    double theta0 = angleOnEllipse(ellipse, x_0);
    double theta1 = angleOnEllipse(ellipse, x_1);
    // Renormalization in case of incoherent angles
    if (theta1 <= theta0) {
        theta1 += 2 * M_PI;
    }
    double result = computeAngularArea(ellipse, theta1) - computeAngularArea(ellipse, theta0);  // only angular part
    // decide whether triangle should be added or removed
    if (theta1 - theta0 < M_PI) {
        result -= geomTools::area::triangle<2>({ 0, 0 }, x_0, x_1);
    } else {
        result += geomTools::area::triangle<2>({ 0, 0 }, x_0, x_1);
    }
    //
    return result;
}

double ellipseAux::computeAngularArea(const Ellipse<2>& ellipse, double theta) {
    //! see https://www.geometrictools.com/Documentation/AreaIntersectingEllipses.pdf
    double a = ellipse.getAlphas()[0];
    double b = ellipse.getAlphas()[1];
    return 0.5 * a * b * (theta -
        atan((b - a) * sin(2 * theta)
            / ((b + a) + (b - a) * cos(2 * theta))));
}

}  // namespace  merope

