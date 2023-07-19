//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef ALGOPACKING_SRC_GEOMETRY_CYLINDER_HXX_
#define ALGOPACKING_SRC_GEOMETRY_CYLINDER_HXX_

#include "../StdHeaders.hxx"

#include "../Geometry/GeomTypes.hxx"
#include "../Geometry/GeomTools_1.hxx"

namespace sac_de_billes {
using namespace std;

//! Class implementing a cylinder
template<unsigned short DIM>
struct Cylinder {
    //! constructor
    Cylinder(const Segment<DIM>& axis_, double radius_) {
        axis = axis_;
        radius = radius_;
    }
    static_assert(DIM == 3);
    // volume of a cylinder
    double volume() { return m_PI * radius * radius * geomTools::norme<DIM>(axis); }
    //! is this point inside?
    bool isInside(const Point<DIM>& point) const;

    //
    Segment<DIM> axis;
    double radius;
    //
    void print(std::ostream& f) const;
};

} // namespace sac_de_billes

#include "../Geometry/Cylinder.ixx"

#endif /* ALGOPACKING_SRC_GEOMETRY_CYLINDER_HXX_ */
