//! Copyright : see license.txt
//!
//! \brief

#ifndef RENORMPOINT_HXX
#define RENORMPOINT_HXX

#include "Point.hxx"
#include "GeomTools_1.hxx"

namespace sac_de_billes {

//! @brief : class implementing a point that is of norm 1
template<unsigned short DIM>
class RenormPoint : private Point<DIM> {
public:
    RenormPoint() : Point<DIM>(create_array<DIM>(0.)) { (*this)(0) = 1.; }
    RenormPoint(const Point<DIM>& pt) : Point<DIM>(pt) { geomTools::renormalize<DIM>(*this); }
    RenormPoint(Point<DIM>&& pt) : Point<DIM>(pt) { geomTools::renormalize<DIM>(*this); }
    inline double operator[](size_t i) const { return Point<DIM>::operator[](i); }
    // nonconst operator[] is deleted for it should be private, but can't.
    // Indeed, it would prevent non const RenormPoint use the const operator.
    const Point<DIM>& getPoint() const { return *this; }
    // friend : lineartransform
    friend double linearTransform::normal<DIM>(RenormPoint<DIM>&, const Point<DIM>&);
private:
    inline double& operator()(size_t i) { return Point<DIM>::operator[](i); }
};

} // namespace sac_de_billes

#endif // RENORMPOINT_HXX
