//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GEOMETRY_ASPRATIO_IXX_
#define GEOMETRY_ASPRATIO_IXX_


#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
inline void WithAspratio<DIM>::setAspRatio(array<double, DIM> aspratio_) {
    auxi_aspratio::renormalizeAspRatio <DIM>(aspratio_);
    aspratio = aspratio_;
    for (size_t i = 0; i < DIM; i++) {
        if (aspratio[i] <= 0) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Nonpositive aspect ratio does not make sense!");
        }
        inverse_aspratio[i] = 1. / aspratio[i];
    }
}

template<unsigned short DIM>
inline void auxi_aspratio::renormalizeAspRatio(Point<DIM>& aspratio) {
    double averaged_geom_aspratio = 1;
    for (size_t i = 0; i < DIM; i++) {
        averaged_geom_aspratio *= aspratio[i];
    }
    averaged_geom_aspratio = pow(averaged_geom_aspratio, 1. / DIM);
    for (size_t i = 0; i < DIM; i++) {
        aspratio[i] /= averaged_geom_aspratio;
    }
}

} // namespace merope

#endif /* GEOMETRY_ASPRATIO_IXX_ */
