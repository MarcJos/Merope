//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../Geometry/Point.hxx"
#include "BasicGeometricOperations.hxx"

namespace sac_de_billes {

template<unsigned short DIM1, unsigned short DIM2>
Point<DIM1>  geomTools::linearProjection(Point<DIM2> oldPoint) {
    static_assert(DIM1 <= DIM2);
    if constexpr (DIM1 == 2 and DIM2 == 3) {
        return Point<DIM1>{oldPoint[0], oldPoint[1]};
    } else {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

inline void geomTools::projection_periodic_1D(double& x, double L) {
    while (x < 0) {
        x += L;
    }
    while (x > L) {
        x -= L;
    }
}

inline void geomTools::projection_periodic_1D_centered(double& x,
    const double& L) {
    while (x < -0.5 * L) {
        x += L;
    }
    while (x > 0.5 * L) {
        x -= L;
    }
}

template<unsigned short DIM>
inline Point<DIM> geomTools::projection(const Segment<DIM>& line,
    const Point<DIM>& z) {
    Point<DIM> xy = line[1] - line[0];
    double segNormCarre = geomTools::normeCarre<DIM>(xy);
    return line[0] + (1. / segNormCarre) * geomTools::prodScal<DIM>(z - line[0], xy) * xy;
}

inline void geomTools::bounce_1D(double& x, const double& L, const double& r) {
    if (x - r < 0) {
        x += 2 * (r - x);
    } else if (x + r > L) {
        x -= 2 * (x + r - L);
    }
}

template<unsigned short DIM, class C1, class C2>
double geomTools::prodScal(const C1& v1, const C2& v2) {
    double res = 0.;
    for (size_t i = 0; i < DIM; i++) {
        res += v1[i] * v2[i];
    }
    return res;
}


template<unsigned short DIM, class C1, class C2>
Point<DIM> geomTools::prodVec(const C1& v1, const C2& v2) {
    if constexpr (DIM == 2) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Incorrect dimension");
    }
    return Point<DIM>{
        v1[1] * v2[2] - v1[2] * v2[1],
            v1[2] * v2[0] - v1[0] * v2[2],
            v1[0] * v2[1] - v1[1] * v2[0]
    };
}

template<unsigned short DIM, class C1, class C2>
inline Point<DIM> geomTools::odot(const C1& v1, const C2& v2) {
    Point<DIM> result;
    for (size_t i = 0; i < DIM; i++) {
        result[i] = v1[i] * v2[i];
    }
    return result;
}

template<unsigned short DIM, class T, typename std::enable_if_t<is_Point<T, DIM>, bool>>
inline double geomTools::normeCarre(const T& vec) {
    return prodScal <DIM>(vec, vec);
}

template<unsigned short DIM, class T, typename std::enable_if_t<is_Point<T, DIM>, bool>>
inline double geomTools::norme(const T& vec) {
    if constexpr (DIM == 1) {
        return abs(vec[0]);
    } else if constexpr (DIM == 2) {
        return hypot(vec[0], vec[1]);
    } else if constexpr (DIM == 3) {
        return hypot(vec[0], vec[1], vec[2]);
    } else {
        return sqrt(normeCarre<DIM>(vec));
    }
}

template<unsigned short DIM, class T, typename std::enable_if_t<is_Point<T, DIM>, bool>>
inline double geomTools::distanceCarre(const T& v1, const T& v2) {
    double res = 0.;
    for (size_t i = 0; i < DIM; i++) {
        res += auxi_function::puissance<2>(v1[i] - v2[i]);
    }
    return res;
}

template<unsigned short DIM, class T, typename std::enable_if_t<is_Point<T, DIM>, bool>>
inline bool geomTools::areEqual(const T& v1, const T& v2, double epsilon) {
    return distanceCarre<DIM>(v1, v2) <= epsilon * (std::numeric_limits<double>::epsilon()
        + geomTools::normeCarre<DIM>(v1) + geomTools::normeCarre<DIM>(v2));
}

template<unsigned short DIM>
inline double geomTools::renormalize(Point<DIM>& v) {
    double norm = sqrt(normeCarre <DIM>(v)) + numeric_limits<double>::min();
    for (size_t i = 0; i < DIM; i++) {
        v[i] /= norm;
    }
    return norm;
}

}  // namespace sac_de_billes


