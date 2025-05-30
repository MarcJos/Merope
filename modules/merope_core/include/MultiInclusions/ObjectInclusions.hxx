//! Copyright : see license.txt
//!
//! \briefClass implementing elliptic inclusions
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"

#include "../MesoStructure/InsideTorus.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"


namespace merope {

template<unsigned short DIM, class OBJECT>
class ObjectInclusions : public InsideTorus<DIM> {
public:
    //! @brief : constructor
    ObjectInclusions() : InsideTorus<DIM>(), objectVector{}{}
    //! @brief : destructor
    virtual ~ObjectInclusions() {}
    ObjectInclusions(array<double, DIM> L_, const vector<OBJECT>& objectVector_) : InsideTorus<DIM>(L_) { this->setInclusions(objectVector_); }
    //! @brief : get all the cells
    const vector<OBJECT>& getMicroInclusions() const { return objectVector; }
    //! @brief : get all the cells
    vector<OBJECT>& getMicroInclusions() { return objectVector; }
    //! @brief : adds another spheroPolyhedron
    void addInclusion(OBJECT objectInclusion) { this->objectVector.emplace_back(objectInclusion); }
    //! @brief : set spheroPolyhedrons
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

template<unsigned short DIM>
using PolyInclusions = ObjectInclusions<DIM, smallShape::ConvexPolyhedronInc<DIM>>;

}  // namespace merope
