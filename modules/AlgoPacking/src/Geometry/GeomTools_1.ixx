//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GEOMTOOLS1_IXX_
#define GEOMTOOLS1_IXX_

namespace sac_de_billes {

template<class VEC_SOLID>
inline double geomTools::volume_all(const VEC_SOLID& s) {
    double res = 0.;
    for (const auto& sph : s) {
        res += geomTools::volume(sph);
    }
    return res;
}

template<unsigned short DIM, class SOLID>
inline bool geomTools::isInside_Intersection(const vector<SOLID>& solids,
    const Point<DIM>& x) {
    for (const auto& solid : solids) {
        if (not isInside<DIM>(solid, x)) {
            return false;
        }
    }
    return true;
}

template<unsigned short DIM, class SOLID>
inline bool sac_de_billes::geomTools::isInside_Union(const vector<SOLID>& solids,
    const Point<DIM>& x) {
    for (const auto& solid : solids) {
        if (isInside<DIM>(solid, x)) {
            return true;
        }
    }
    return false;
}

template<unsigned short DIM, class SOLID>
inline double geomTools::distanceTo(const SOLID& solid,
    const Point<DIM>& vector_from_center_to_point) {
    // CASE 1 : halspaces
    if constexpr (is_same<SOLID, std::vector<HalfSpace<DIM>>>::value) {
        double distance = -0.5 * numeric_limits<double>::max();;
        double currentDistance = 0;
        for (const auto& hspace : solid) {
            currentDistance = hspace.distanceTo(vector_from_center_to_point);
            distance = max(currentDistance, distance);
        }
        return distance;
    }
    // CASE 2 : polyhedron
    else if constexpr (is_same<SOLID, ConvexPolyhedron<DIM>>::value) {
        return distanceTo<DIM>(solid.faces, vector_from_center_to_point);
    }
    // CASE 3 : sphere
    else if constexpr (is_same<SOLID, Sphere<DIM>>::value) {
        return geomTools::norme<DIM>(vector_from_center_to_point) - solid.radius;
    }
    // CASE 4 : ellipse
    else if constexpr (is_same<SOLID, Ellipse<DIM>>::value) {
        // should be modified by Leo
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected");
    }
    // ELSE : not programmed
    else {
        throw runtime_error("Not programemd");
    }
}

/// RENORMALIZE

template<unsigned short DIM>
inline void linearTransform::point(Point<DIM>& vec,
    const Point<DIM>& linTransform) {
    for (auto i = 0; i < DIM; i++) {
        vec[i] *= linTransform[i];
    }
}

template<unsigned short DIM, class C>
inline double linearTransform::normal(C& normal,
    const Point<DIM>& linTransform) {
    for (auto i = 0; i < DIM; i++) {
        normal(i) /= linTransform[i];
    }
    return geomTools::renormalize <DIM>(normal);
}

template<unsigned short DIM>
inline double linearTransform::determinant(const Point<DIM>& linTransform) {
    double res = 1;
    for (size_t i = 0; i < DIM; i++) {
        res *= linTransform[i];
    }
    return res;
}

template<unsigned short DIM>
inline void linearTransform::proceed(HalfSpace<DIM>& halfspace,
    const Point<DIM>& linTransform) {
    double norm = linearTransform::normal <DIM>(halfspace.vec_force_definition(), linTransform);
    halfspace.c() /= norm;
}

template<unsigned short DIM>
inline void linearTransform::proceed(Sphere<DIM>& sphere,
    const Point<DIM>& linTransform) {
    linearTransform::point<DIM>(sphere.center, linTransform);
    cerr << __PRETTY_FUNCTION__ << endl;
    throw runtime_error("Incomplete");
}

template<unsigned short DIM>
inline void  linearTransform::proceed(
    ConvexPolyhedron<DIM>& convexPolyhedron,
    const Point<DIM>& linTransform) {
    linearTransform::point <DIM>(convexPolyhedron.center, linTransform);
    for (auto& f : convexPolyhedron.faces) {
        linearTransform::proceed <DIM>(f, linTransform);
    }
}

template<unsigned short DIM>
inline void linearTransform::proceed(Ellipse<DIM>& ellipse,
    const Point<DIM>& linTransform) {
    cerr << __PRETTY_FUNCTION__ << endl;
    throw runtime_error("Incomplete");
}


} // namespace sac_de_billes

#endif /* GEOMTOOLS_IXX_ */
