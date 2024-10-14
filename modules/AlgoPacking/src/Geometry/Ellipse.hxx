//! Copyright : see license.txt
//!
//! \briefClass for defining ellipses
//
#pragma once

#include "../StdHeaders.hxx"

#include "../Geometry/GeomTypes.hxx"
#include "../Geometry/GeomTools_1.hxx"

namespace sac_de_billes {
using namespace std;

template<unsigned short DIM>
class Ellipse {
    //! x inside ellipse if \sum |x \cdot a_i|^2 < 1
public:
    //! constructor
    Ellipse();
    //! constructor
    explicit Ellipse(const Sphere<DIM> sphere);
    //! constructor
    Ellipse(const Point<DIM>& center_, const array<Point<DIM>, DIM>& axes_);
    //! computes the volume of an ellipse
    double volume() const;
    //! is this point inside?
    bool isInside(const Point<DIM>& point) const;
    //! getter
    const array<double, DIM>& getAlphas() const { return alphas; }
    //! print
    void print(ostream& os) const;

private:
    //! compute the quantities alphas
    void compute_alphas();

public:
    //! phase
    PhaseType phase;
    //! center
    Point<DIM> center;
    //! normed axes of the ellipse, a_i
    array<Point<DIM>, DIM> axes;

private:
    //! return 1/\alpha_i = |a_i|^2
    array<double, DIM> alphas;
};

namespace ellipseAux {
template<unsigned short DIM>
array<Point<DIM>, DIM> defaultAxes();
//! @return the area delimitated by the chord on the ellipse of extremities x_0 and x_1 (assumed to be on the ellipse) [in the direct sense]
//! @param ellipse : ellipse
//! @param x_0 : 1st point delimitating the chord
//! @param x_1  : 2nd point delimitating the chord
double computeChordArea(const Ellipse<2>& ellipse, Point<2> x_0, Point<2> x_1);
//! @return : a polar representation of a point on the ellipse
//! @param ellipse : ellipse
//! @param pointOnEllipse : point assumed to be on the ellipse
double angleOnEllipse(const Ellipse<2>& ellipse, Point<2> pointOnEllipse);
//! @return the area delimitated by the angular portion of the ellipse delimitated by the first axis and of angle theta
//! @param ellipse : ellipse
//! @param theta : angle of the angular portion, starting from the first axis
double computeAngularArea(const Ellipse<2>& ellipse, double theta);
}  // namespace  ellipseAux

}  // namespace sac_de_billes

#include "../Geometry/Ellipse.ixx"

