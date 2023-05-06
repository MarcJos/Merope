//! Copyright : see license.txt
//!
//! \brief
//
#ifndef GLOBALSHAPE_IXX_
#define GLOBALSHAPE_IXX_

///--------------------------------------
/// BigShape
///--------------------------------------

namespace sac_de_billes {
template<unsigned short DIM>
inline double AmbiantSpace::BigShape<DIM>::normeCarre(
        const Point<DIM> &vec) const {
    double result = 0.;
    for (size_t i = 0; i < DIM; i++) {
        result += vec[i] * vec[i];
    }
    return result;
}

template<unsigned short DIM>
inline double AmbiantSpace::BigShape<DIM>::norme(const Point<DIM> &vec) const {
    if constexpr (DIM == 1) {
        return abs(vec[0]);
    } else if constexpr (DIM == 2) {
        return hypot(vec[0], vec[1]);
    } else if constexpr (DIM == 3) {
        return hypot(vec[0], vec[1], vec[2]);
    } else {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM>
bool AmbiantSpace::BigShape<DIM>::areSphereIntersected(const Sphere<DIM> &sph1,
        const Sphere<DIM> &sph2) const {
    return distanceCarre(sph1.center, sph2.center)
            < (sph1.radius + sph2.radius) * (sph1.radius + sph2.radius);
}

///--------------------------------------
/// Cube
///--------------------------------------

template<unsigned short DIM>
AmbiantSpace::Cube<DIM>::Cube(Point<DIM> L_) :
        AmbiantSpace::BigShape<DIM>(L_) {
}

template<unsigned short DIM>
double AmbiantSpace::Cube<DIM>::volume() const {
    if constexpr (DIM == 1) {
        return this->L[0];
    } else if constexpr (DIM == 2) {
        return this->L[0] * this->L[1];
    } else if constexpr (DIM == 3) {
        return this->L[0] * this->L[1] * this->L[2];
    } else {
        throw invalid_argument(
                "Cube<DIM>::volume() La dimension est 1, 2, ou 3");
    }
}

template<unsigned short DIM>
bool AmbiantSpace::Cube<DIM>::isInside(const Point<DIM> &point,
        const double &minRadius) const {
    for (unsigned short i = 0; i < DIM; i++) {
        if (point[i] + minRadius + this->boundaryExclusionDistance > this->L[i]
                or point[i] - minRadius - this->boundaryExclusionDistance < 0) {
            return false;
        }
    }
    return true;
}

template<unsigned short DIM>
inline bool AmbiantSpace::Cube<DIM>::isInside(const Sphere<DIM> &sph) const {
    return isInside(sph.center, sph.radius);
}

template<unsigned short DIM>
double AmbiantSpace::Cube<DIM>::distanceCarre(const Point<DIM> &x1,
        const Point<DIM> &x2) const {
    double result = 0.;
    double dx = 0.;
    for (size_t i = 0; i < DIM; i++) {
        dx = x1[i] - x2[i];
        result += dx * dx;
    }
    return result;
}

template<unsigned short DIM>
Point<DIM> AmbiantSpace::Cube<DIM>::geomVector(const Point<DIM> &p1,
        const Point<DIM> &p2) const {
    Point <DIM> result { };
    for (size_t i = 0; i < DIM; i++) {
        result[i] = p2[i] - p1[i];
    }
    return result;
}

template<unsigned short DIM>
inline void AmbiantSpace::Cube<DIM>::projection(Point<DIM> &point) const {
    for (size_t i = 0; i < DIM; i++) {
        if (point[i] < 0.) {
            point[i] = 0.;
        } else if (point[i] > this->L[i]) {
            point[i] = this->L[i];
        }
    }
}

template<unsigned short DIM>
inline void AmbiantSpace::Cube<DIM>::bounce(Sphere<DIM> &sph) const {
    for (size_t i = 0; i < DIM; i++) {
        geomTools::bounce_1D(sph.center[i], this->L[i],
                sph.radius + this->boundaryExclusionDistance);
    }
}

template<unsigned short DIM>
inline double AmbiantSpace::Cube<DIM>::faceArea(int direction) {
    if (direction < 0 or direction > DIM){
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected");
    }
    double result = 1;
    for(int i = 0; i < DIM; i++){
        if (i != direction) result *= this->L[i];
    }
    return result;
}
///--------------------------------------
/// Tore
///--------------------------------------

template<unsigned short DIM>
AmbiantSpace::Tore<DIM>::Tore(Point<DIM> L_) :
        Cube<DIM>(L_) {
}

template<unsigned short DIM>
double AmbiantSpace::Tore<DIM>::distanceCarre(const Point<DIM> &x1,
        const Point<DIM> &x2) const {
    return this->normeCarre(geomVector(x1, x2));
}

template<unsigned short DIM>
bool AmbiantSpace::Tore<DIM>::isInside(const Point<DIM> &point,
        const double&) const {
    for (unsigned short i = 0; i < DIM; i++) {
        if (point[i] > this->L[i] or point[i] < 0) {
            return false;
        }
    }
    return true;
}

template<unsigned short DIM>
inline bool AmbiantSpace::Tore<DIM>::isInside(const Sphere<DIM>&) const {
    return true;
}

template<unsigned short DIM>
inline void AmbiantSpace::Tore<DIM>::projection(Point<DIM> &point) const {
    for (size_t i = 0; i < DIM; i++) {
        geomTools::projection_periodic_1D(point[i], this->L[i]);
    }
}

template<unsigned short DIM>
Point<DIM> AmbiantSpace::Tore<DIM>::geomVector(const Point<DIM> &p1,
        const Point<DIM> &p2) const {
    Point <DIM> result { };
    for (size_t i = 0; i < DIM; i++) {
        result[i] = p2[i] - p1[i];
        geomTools::projection_periodic_1D_centered(result[i], this->L[i]);
    }
    return result;
}

template<unsigned short DIM>
inline void AmbiantSpace::Tore<DIM>::bounce(Sphere<DIM> &sph) const {
    projection(sph.center);
}

///--------------------------------------
/// BigSphere
///--------------------------------------
template<unsigned short DIM>
inline AmbiantSpace::BigSphere<DIM>::BigSphere(Point<DIM> L_) :
        AmbiantSpace::Cube<DIM>(L_) {
    radius = 0.5 * (*min_element(L_.begin(), L_.end()));
    for (size_t i = 0; i < DIM; i++) {
        center[i] = radius;
    }
}

template<unsigned short DIM>
inline bool AmbiantSpace::BigSphere<DIM>::isInside(const Point<DIM> &point,
        const double &minRadius) const {
    // \fixme if minRadius < Radius, problem.
    return this->distanceCarre(center, point) < auxi_function::puissance < 2
            > (radius - minRadius - this->boundaryExclusionDistance);
}

template<unsigned short DIM>
inline double AmbiantSpace::BigSphere<DIM>::volume() const {
    return sphereTools::volumeSphere <DIM> (radius);
}

template<unsigned short DIM>
inline bool AmbiantSpace::BigSphere<DIM>::isInside(
        const Sphere<DIM> &sph) const {
    return isInside(sph.center, sph.radius);
}

template<unsigned short DIM>
inline void AmbiantSpace::BigSphere<DIM>::bounce(Sphere<DIM> &sph) const {
    double R = sqrt(this->distanceCarre(sph.center, center));
    double delta = R + sph.radius + this->boundaryExclusionDistance - radius;
    if (delta > 0) {
        for (size_t i = 0; i < DIM; i++) {
            sph.center[i] -= center[i]; // computes the sphere wrt the basis with the center of the BigSphere as origin
            sph.center[i] *= (1 - 2 * delta / R); 	// bounce on the boundary
            sph.center[i] += center[i]; // puts the sphere back in the original basis
        }
    }
}

///--------------------------------------
/// BigCylinder
///--------------------------------------
template<unsigned short DIM>
inline AmbiantSpace::BigCylinder<DIM>::BigCylinder(Point<DIM> L_) :
        AmbiantSpace::Cube<DIM>(L_), circle(
                AmbiantSpace::BigSphere < 2
                        > (array<double, 2> { L_[0], L_[1] })) {
    circle.setBoundaryExclusionDistance(this->boundaryExclusionDistance);
    if (DIM != 3) {
        throw invalid_argument(__PRETTY_FUNCTION__);
    }
    height = L_[2];
}

template<unsigned short DIM>
inline bool AmbiantSpace::BigCylinder<DIM>::isInside(const Point<DIM> &point,
        const double &minRadius) const {
    // check the height
    if (point[2] + minRadius + this->boundaryExclusionDistance > height
            or point[2] - (minRadius + this->boundaryExclusionDistance) < 0) {
        return false;
    }
    // check the circle
    Point < 2 > point_2dim = { point[0], point[1] };
    return circle.isInside(point_2dim, minRadius);
}

template<unsigned short DIM>
inline double AmbiantSpace::BigCylinder<DIM>::volume() const {
    return circle.volume() * height;
}

template<unsigned short DIM>
inline bool AmbiantSpace::BigCylinder<DIM>::isInside(
        const Sphere<DIM> &sph) const {
    return isInside(sph.center, sph.radius);
}

template<unsigned short DIM>
inline void AmbiantSpace::BigCylinder<DIM>::bounce(Sphere<DIM> &sph) const {
    // radial part
    Sphere <2> circleSph(Point<2> { sph.center[0], sph.center[1] }, sph.radius, sph.phase);
    circle.bounce(circleSph);
    sph.center[0] = circleSph.center[0];
    sph.center[1] = circleSph.center[1];
    // vertical part
    geomTools::bounce_1D(sph.center[2], this->L[2],
            sph.radius + this->boundaryExclusionDistance);
}

} // namespace sac_de_billes

#endif /* GLOBALSHAPE_IXX_ */
