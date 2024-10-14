//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "GeomTools_1.hxx"

namespace sac_de_billes {
namespace geomTools {
namespace area {

template<unsigned short DIM>
double triangle(Point<DIM> pt0, Point<DIM> pt1, Point<DIM> pt2) {
    double a = geomTools::norme<DIM>(pt1 - pt0);
    double b = geomTools::norme<DIM>(pt2 - pt1);
    double c = geomTools::norme<DIM>(pt0 - pt2);
    double p = 0.5 * (a + b + c);
    return sqrt(p * (p - a) * (p - b) * (p - c));
}

}  // namespace area
}  // namespace geomTools
}  // namespace  sac_de_billes



