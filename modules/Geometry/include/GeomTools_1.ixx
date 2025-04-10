//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace merope {

template<class VEC_SOLID>
double geomTools::volume_all(const VEC_SOLID& s) {
    double res = 0.;
    for (const auto& sph : s) {
        res += geomTools::volume(sph);
    }
    return res;
}

template<unsigned short DIM, class SOLID>
bool geomTools::isInside_Intersection(const vector<SOLID>& solids,
    const Point<DIM>& x) {
    for (const auto& solid : solids) {
        if (not isInside<DIM>(solid, x)) {
            return false;
        }
    }
    return true;
}

template<unsigned short DIM, class SOLID>
bool merope::geomTools::isInside_Union(const vector<SOLID>& solids,
    const Point<DIM>& x) {
    for (const auto& solid : solids) {
        if (isInside<DIM>(solid, x)) {
            return true;
        }
    }
    return false;
}

template<unsigned short DIM, class SOLID>
double geomTools::distanceTo(const SOLID& solid,
    const Point<DIM>& vector_from_center_to_point) {
    // CASE 1 : halspaces
    if constexpr (is_same_v<SOLID, std::vector<HalfSpace<DIM>>>) {
        double distance = -0.5 * numeric_limits<double>::max();;
        double currentDistance = 0;
        for (const auto& hspace : solid) {
            currentDistance = hspace.distanceTo(vector_from_center_to_point);
            distance = max(currentDistance, distance);
        }
        return distance;
    }
    // CASE 2 : polyhedron
    else if constexpr (is_same_v<SOLID, ConvexPolyhedron<DIM>>) {
        return distanceTo<DIM>(solid.faces, vector_from_center_to_point);
    }
    // CASE 3 : sphere
    else if constexpr (is_same_v<SOLID, Sphere<DIM>>) {
        return geomTools::norme<DIM>(vector_from_center_to_point) - solid.radius;
    }
    // CASE 4 : ellipse
    else if constexpr (is_same_v<SOLID, Ellipse<DIM>>) {
        Merope_error_not_done();
    }
    // CASE 5 : SpheroPolyhedron
    else if constexpr (is_same_v<SOLID, SpheroPolyhedron<DIM>>) {
        return solid.distanceTo(vector_from_center_to_point);
    }
    // ELSE : not programmed
    else {
        Merope_error_not_done();
    }
}

/// RENORMALIZE

template<unsigned short DIM>
void linearTransform::point(Point<DIM>& vec,
    const Point<DIM>& linTransform) {
    for (auto i = 0; i < DIM; i++) {
        vec[i] *= linTransform[i];
    }
}

template<unsigned short DIM, class C>
double linearTransform::normal(C& normal,
    const Point<DIM>& linTransform) {
    for (auto i = 0; i < DIM; i++) {
        normal(i) /= linTransform[i];
    }
    return geomTools::renormalize <DIM>(normal);
}

template<unsigned short DIM>
double linearTransform::determinant(const Point<DIM>& linTransform) {
    double res = 1;
    for (size_t i = 0; i < DIM; i++) {
        res *= linTransform[i];
    }
    return res;
}

template<unsigned short DIM>
void linearTransform::proceed(HalfSpace<DIM>& halfspace,
    const Point<DIM>& linTransform) {
    double norm = linearTransform::normal <DIM>(halfspace.vec_force_definition(), linTransform);
    halfspace.c() /= norm;
}

template<unsigned short DIM>
void linearTransform::proceed(Sphere<DIM>& sphere,
    const Point<DIM>& linTransform) {
    linearTransform::point<DIM>(sphere.center, linTransform);
    cerr << __PRETTY_FUNCTION__ << endl;
    throw runtime_error("Incomplete");
}

template<unsigned short DIM>
void  linearTransform::proceed(
    ConvexPolyhedron<DIM>& convexPolyhedron,
    const Point<DIM>& linTransform) {
    linearTransform::point <DIM>(convexPolyhedron.center, linTransform);
    for (auto& f : convexPolyhedron.faces) {
        linearTransform::proceed <DIM>(f, linTransform);
    }
}

template<unsigned short DIM>
void linearTransform::proceed(Ellipse<DIM>& ellipse,
    const Point<DIM>& linTransform) {
    cerr << __PRETTY_FUNCTION__ << endl;
    throw runtime_error("Incomplete");
}


}  // namespace merope


