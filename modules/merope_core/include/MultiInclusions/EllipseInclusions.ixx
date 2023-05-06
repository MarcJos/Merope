//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MESOSTRUCTURE_ELLIPSEINCLUSIONS_IXX_
#define MESOSTRUCTURE_ELLIPSEINCLUSIONS_IXX_

#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
EllipseInclusions<DIM>::EllipseInclusions(array<double, DIM> L_, const vector<Ellipse<DIM>>& ellipses_):
    InsideTorus<DIM>(L_), ellipses(ellipses_) {}

} // namespace merope


#endif /* MESOSTRUCTURE_ELLIPSEINCLUSIONS_IXX_ */
