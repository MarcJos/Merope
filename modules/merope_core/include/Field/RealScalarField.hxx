//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/GeomTypes.hxx"

namespace merope {
namespace realScalarField {

template<unsigned short DIM>
//! just for plotting a real field
class Field {
public:
    //! constructor
    Field(std::function<double(Point<DIM>)> fieldFunction_) : fieldFunction{ fieldFunction_ } {};
    //! stores the fieldFunction
    std::function<double(Point<DIM>)> fieldFunction;
};

}  // namespace realScalarField
}  // namespace merope



