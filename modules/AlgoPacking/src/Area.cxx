//! Copyright : see license.txt
//!
//! \brief

#include "Geometry/Area.hxx"
#include "Geometry/BasicGeometricOperations.hxx"


double sac_de_billes::geomTools::area::polygon(vector<Point<2>> vertices) {
    size_t n = vertices.size();
    double result = geomTools::determinant(vertices[n - 1], vertices[0]);
    for (size_t i = 0; i + 1 < n; i++) {
        result += geomTools::determinant(vertices[i], vertices[i + 1]);
    }
    return 0.5 * result;
}


