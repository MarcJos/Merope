//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MEROPE_CORE_SRC_VOXELLATION_REALSCALARFIELD_HXX_
#define MEROPE_CORE_SRC_VOXELLATION_REALSCALARFIELD_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/Geometry/GeomTypes.hxx"



#include "../MeropeNamespace.hxx"


namespace merope {
namespace realScalarField {

template<unsigned short DIM>
//! just for plotting a real field
class Field {
public:
    //! constructor
    Field(std::function<double(Point<DIM>)> fieldFunction_): fieldFunction{ fieldFunction_ } {};
    //! stores the fieldFunction
    std::function<double(Point<DIM>)> fieldFunction;
};

} // namespace realScalarField
} // namespace merope


#endif /* MEROPE_CORE_SRC_VOXELLATION_REALSCALARFIELD_HXX_ */
