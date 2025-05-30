//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "GeomTypes.hxx"
#include "GeomTools_1.hxx"

namespace merope {

//! Class implementing a cylinder
template<unsigned short DIM, typename = std::enable_if_t<DIM == 3, bool>>
struct Cylinder {
    //! constructor
    Cylinder(const Segment<DIM>& axis_, double radius_, PhaseType phase_ = 0) :
        axis{ axis_ }, radius{ radius_ }, phase{ phase_ } {}
    // volume of a cylinder
    double volume() const { return m_PI * radius * radius * geomTools::norme<DIM>(axis[1] - axis[0]); }
    //! is this point inside?
    bool isInside(const Point<DIM>& point) const;
    Point<DIM> center() const { return 0.5 * (axis[0] + axis[1]); }
    //
    Segment<DIM> axis;
    double radius;
    PhaseType phase;
    //
    void print(std::ostream& f) const;
};

}  // namespace merope

#include "Cylinder.ixx"


