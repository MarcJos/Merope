//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Geometry/GeomTools.hxx"

#include "../MeropeNamespace.hxx"


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


