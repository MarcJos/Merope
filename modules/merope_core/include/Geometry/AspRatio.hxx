//! Copyright : see license.txt
//!
//! \briefImplement an aspect ratio shared by the whole structure
//
//

#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
class WithAspratio {
public:
    //! \param aspratio_ (optional) : aspect ratio
    void setAspRatio(array<double, DIM> aspratio_);
protected:
    WithAspratio() {
        setAspRatio(create_array<DIM>(1.));
    }
    //! aspect ratio (renormalized)
    array<double, DIM> aspratio;
    //! inverse of aspect ratio (renormalized)
    array<double, DIM> inverse_aspratio;
};


namespace auxi_aspratio {
template<unsigned short DIM>
void renormalizeAspRatio(Point<DIM>&);
}  // auxi_aspratio

}  // namespace merope

#include "../Geometry/AspRatio.ixx"


