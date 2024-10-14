//! Copyright : see license.txt
//!
//! \briefClass implementing elliptic inclusions
//
#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../MesoStructure/InsideTorus.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM, class OBJECT>
class ObjectInclusions final : public InsideTorus<DIM> {
public:
    //! constructor
    ObjectInclusions() : InsideTorus<DIM>(), objectVector{}{}
    ObjectInclusions(array<double, DIM> L_, const vector<OBJECT>& objectVector_) : InsideTorus<DIM>(L_) { this->setInclusions(objectVector_); }
    //! get all the cells
    const vector<OBJECT>& getMicroInclusions() const { return objectVector; }
    //! adds another spheroPolyhedron
    void addInclusion(OBJECT objectInclusion) { this->objectVector.emplace_back(objectInclusion); }
    //! set spheroPolyhedrons
    void setInclusions(const vector<OBJECT>& objectVector_) { this->objectVector = objectVector_; }

protected:
    //! spheroPolyhedrons
    vector<OBJECT> objectVector;
};

template<unsigned short DIM>
using SpheroPolyhedronInclusions = ObjectInclusions<DIM, smallShape::SpheroPolyhedronInc<DIM>>;

template<unsigned short DIM>
using EllipseInclusions = ObjectInclusions<DIM, smallShape::EllipseInc<DIM>>;

template<unsigned short DIM>
using CylinderInclusions = ObjectInclusions<DIM, smallShape::CylinderInc<DIM>>;

}  // namespace merope
