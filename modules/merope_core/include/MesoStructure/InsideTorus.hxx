//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"
#include "../../../Geometry/include/GeomTools.hxx"


namespace merope {

template<unsigned short DIM>
class InsideTorus {
private:
    static constexpr array<double, DIM> L_DEFAULT = create_array<DIM>(1.);
public:
    //! constructor
    InsideTorus(array<double, DIM> L) :
        tore{ AmbiantSpace::Tore<DIM>(L) } {}
    //! constructor
    InsideTorus() : InsideTorus<DIM>(L_DEFAULT) {}
    //! destructor
    virtual ~InsideTorus() {}

    //! get the length
    const array<double, DIM>& getL() const {
        return tore.L;
    }
    //! Set length
    virtual void setLength(array<double, DIM> length) {
        tore = AmbiantSpace::Tore<DIM>(length);
    }

    //! Ambiant space
    AmbiantSpace::Tore<DIM> tore;
};

}  // namespace merope


