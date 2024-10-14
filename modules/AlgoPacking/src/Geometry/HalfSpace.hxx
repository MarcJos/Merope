//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../Geometry/GeomTypes.hxx"
#include "../Geometry/GeomConstants.hxx"
#include "../Geometry/RenormPoint.hxx"

namespace sac_de_billes {
using namespace std;

template<unsigned short DIM>
//! x \in \R^DIM such that x \cdot vec < c
//! vec is assumed to be of norm 1
class HalfSpace {
public:
    //! constructor
    HalfSpace(Point<DIM> vec__, double c__);
    //! constructor
    HalfSpace(Point<DIM> outerNormal, Point<DIM> pointPlane);

    //! \return whether a point is inside the half-space
    bool isInside(const Point<DIM>& x) const;
    //! \return the signed distance of a point w.r.t to the surface halfspace (<0 = inside, >0 = outside)
    double distanceTo(const Point<DIM>& x)const;

    //! for showing it
    string toString();

    //! getter
    double& c() { return c_; }
    //! getter
    double c() const { return c_; }
    //! getter
    const Point<DIM>& vec() const { return vec_.getPoint(); }
    //! \warning should preserve the norm = 1
    RenormPoint<DIM>& vec_force_definition() { return vec_; }

private:
    double c_;
    //! vec is forced to be of norm 1
    RenormPoint<DIM> vec_;
};

}  // namespace sac_de_billes

#include "../Geometry/HalfSpace.ixx"


