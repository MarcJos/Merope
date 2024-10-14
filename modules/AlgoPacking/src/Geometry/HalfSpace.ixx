//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../Geometry/GeomTools_1.hxx"

namespace sac_de_billes {

///HalfSpace<DIM>

template<unsigned short DIM>
HalfSpace<DIM>::HalfSpace(Point<DIM> vec__, double c__) :
    c_{ c__ }, vec_{ vec__ } {}

template<unsigned short DIM>
HalfSpace<DIM>::HalfSpace(Point<DIM> outerNormal, Point<DIM> pointPlane) :
    HalfSpace(outerNormal, geomTools::prodScal<DIM>(outerNormal, pointPlane)) {}

template<unsigned short DIM>
bool HalfSpace<DIM>::isInside(const Point<DIM>& x) const {
    return geomTools::prodScal<DIM>(x, vec_) < c_;
}

template<unsigned short DIM>
double HalfSpace<DIM>::distanceTo(const Point<DIM>& x)const {
    return geomTools::prodScal<DIM>(x, vec_) - c_;
}

template<unsigned short DIM>
inline string HalfSpace<DIM>::toString() {
    string result = "HalfSpace \n";
    result += "c = " + to_string(c_) + "\n";
    result += "vec = (";
    for (size_t i = 0; i < DIM; i++) {
        result += to_string(vec_[i]) + "   ";
    }
    result += ")\n";
    return result;
}

}  // namespace sac_de_billes


