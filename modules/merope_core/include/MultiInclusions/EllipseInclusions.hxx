//! Copyright : see license.txt
//!
//! \brief Class implementing elliptic inclusions
//
#ifndef MESOSTRUCTURE_ELLIPSEINCLUSIONS_HXX_
#define MESOSTRUCTURE_ELLIPSEINCLUSIONS_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../MesoStructure/InsideTorus.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
class EllipseInclusions: public InsideTorus<DIM> {
    // Class implementing elliptic inclusions
public:
    //! \param L_ : the torus dimensions
    //! \param seeds_  : a list of seeds parametrizing the tessellation (=centers of tessels + weight + phase Id)
    EllipseInclusions(array<double, DIM> L_, const vector<Ellipse<DIM>>& ellipses_);
    //! get all the cells
    vector<smallShape::EllipseInc<DIM>> getMicroInclusions() const { return ellipses; }
    //! ellipses
    vector<smallShape::EllipseInc<DIM>> ellipses;
};

} // namespace merope


#include "../MultiInclusions/EllipseInclusions.ixx"

#endif /* MESOSTRUCTURE_ELLIPSEINCLUSIONS_HXX_ */
