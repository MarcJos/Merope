//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MESOSTRUCTURE_POLYINCLUSIONS_HXX_
#define MESOSTRUCTURE_POLYINCLUSIONS_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AlgoNames.hxx"
#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Geometry/AspRatio.hxx"
#include "../MesoStructure/InsideTorus.hxx"
#include "../MeropeNamespace.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"


namespace merope {

template<unsigned short DIM>
class PolyInclusions: public InsideTorus<DIM>, public WithAspratio<DIM> {
    // class implementing non-intersecting polyhedric inclusions
public:
    //! constructor
    PolyInclusions(): InsideTorus<DIM>(), WithAspratio<DIM>(), polyInclusions{}{};
    //! destructor
    virtual ~PolyInclusions() {};
    //! get all the cells
    const vector<smallShape::ConvexPolyhedronInc<DIM>>& getMicroInclusions() const { return polyInclusions; }
    //! adds a new polyhedron
    virtual void addInclusion(const smallShape::ConvexPolyhedronInc<DIM>& polyhedron) { this->polyInclusions.push_back(polyhedron); };
    //! set spheroPolyhedrons
    void setInclusions(vector<smallShape::ConvexPolyhedronInc<DIM>> polyInclusions_) { this->polyInclusions = polyInclusions_; }

protected:
    //! tessels
    vector<smallShape::ConvexPolyhedronInc<DIM>> polyInclusions;
};

} // namespace merope


#include "../MultiInclusions/PolyInclusions.ixx"

#endif /* MESOSTRUCTURE_POLYINCLUSIONS_HXX_ */
