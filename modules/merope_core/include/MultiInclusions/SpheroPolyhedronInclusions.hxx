//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MEROPE_CORE_SRC_MULTIINCLUSIONS_SPHEROPOLYHEDRONINCLUSIONS_HXX_
#define MEROPE_CORE_SRC_MULTIINCLUSIONS_SPHEROPOLYHEDRONINCLUSIONS_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"
#include "../../../AlgoPacking/src/AmbiantSpace.hxx"

#include "../Geometry/GeomTools.hxx"
#include "../MesoStructure/InsideTorus.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"

#include "../MeropeNamespace.hxx"

namespace merope {

//! class implementing non-intersecting spheroPolyhedric inclusions
template<unsigned short DIM>
class SpheroPolyhedronInclusions: public InsideTorus<DIM> {
public:
    //! constructor
    SpheroPolyhedronInclusions(): InsideTorus<DIM>() {};
    //! get all the cells
    const vector<smallShape::SpheroPolyhedronInc<DIM>>& getMicroInclusions() const { return spheroPolyInclusions; }
    //! adds another spheroPolyhedron
    void addInclusion(smallShape::SpheroPolyhedronInc<DIM> spheroPolyInclusion) { this->spheroPolyInclusions.emplace_back(spheroPolyInclusion); }
    //! set spheroPolyhedrons
    void setInclusions(vector<smallShape::SpheroPolyhedronInc<DIM>> spheroPolyInclusions_) { this->spheroPolyInclusions = spheroPolyInclusions_; }

protected:
    //! spheroPolyhedrons
    vector<smallShape::SpheroPolyhedronInc<DIM>> spheroPolyInclusions;
};


} // namespace merope

#endif /* MEROPE_CORE_SRC_MULTIINCLUSIONS_SPHEROPOLYHEDRONINCLUSIONS_HXX_ */
