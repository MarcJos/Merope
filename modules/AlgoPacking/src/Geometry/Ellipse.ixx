#include "Ellipse.hxx"
//! Copyright : see license.txt
//!
//! \brief
//
#ifndef ELLIPSE_IXX_
#define ELLIPSE_IXX_

namespace sac_de_billes {

template<unsigned short DIM>
Ellipse<DIM>::Ellipse(): Ellipse(Sphere<DIM>()) {}

template<unsigned short DIM>
Ellipse<DIM>::Ellipse(const Sphere<DIM> sphere) : Ellipse(sphere.center, ellipseAux::defaultAxes<DIM>()) {
    for (auto& ax : this->axes) {
        ax *= 1. / sphere.radius;
    }
    compute_alphas();
}

template<unsigned short DIM>
inline Ellipse<DIM>::Ellipse(const Point<DIM>& center_,
    const array<Point<DIM>, DIM>& axes_): phase{ 0 }, center(center_), axes(axes_) {
    compute_alphas();
}

template<unsigned short DIM>
inline void Ellipse<DIM>::compute_alphas() {
    for (size_t i = 0; i < DIM; i++) {
        alphas[i] = 1. / geomTools::norme<DIM>(axes[i]);
    }
}

template<unsigned short DIM>
inline double Ellipse<DIM>::volume() const {
    if constexpr (DIM == 2) {
        return m_PI * alphas[0] * alphas[1];
    } else if constexpr (DIM == 3) {
        return m_PI * 4. / 3. * alphas[0] * alphas[1] * alphas[2];
    } else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected!");
        //static_assert(false);
    }
}

template<unsigned short DIM>
inline bool Ellipse<DIM>::isInside(const Point<DIM>& point) const {
    double res = 0;
    for (size_t i = 0; i < DIM; i++) {
        res += auxi_function::puissance<2>(geomTools::prodScal<DIM>(point - this->center, axes[i]));
    }
    return res < 1;
}


template<unsigned short DIM>
array<Point<DIM>, DIM> ellipseAux::defaultAxes() {
    array<Point<DIM>, DIM> result = create_array<DIM, Point<DIM>>(create_array < DIM>(0.));
    for (size_t i = 0; i < DIM; i++) {
        result[i][i] = 1;
    }
    return result;
}

} // namespace sac_de_billes

#endif /* ELLIPSE_IXX_ */
