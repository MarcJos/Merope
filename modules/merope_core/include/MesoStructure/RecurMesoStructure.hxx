//! Copyright : see license.txt
//!
//! \brief
#include "../SingleVoxel/RecurStructure.hxx"

#include "../Field/CartesianField.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"

#pragma once

namespace merope {

template<unsigned short DIM, class BasicStruct, class BasicType>
class RecursiveStructure_inTorus : public RecursiveStructure<DIM, BasicStruct, BasicType> {
    static_assert(std::is_base_of<InsideTorus<DIM>, BasicStruct>::value);
public:
    const array<double, DIM>& getL() const {
        return this->getOneBasicStructure().getL();
    }
};


template<unsigned short DIM>
using RecurField = RecursiveStructure_inTorus<DIM, CartesianField<DIM>, double>;

template<unsigned short DIM>
class RecurGeom : public RecursiveStructure_inTorus<DIM, MultiInclusions<DIM>, PhaseType> {
    using BasicStruct = MultiInclusions<DIM>;
    using BasicType = PhaseType;
public:
    //! Default initialization
    RecurGeom();
    //! get all contained phases
    const vector<PhaseType>& getAllPhases() const { return allPhases; }
    //! setter
    void setMultiInclusions(const MultiInclusions<DIM>& mi);
    void setCombination2(const RecurGeom<DIM>& mi1, const RecurGeom<DIM>& mi2,
        std::function<PhaseType(PhaseType, PhaseType)> func,
        vector<PhaseType> allPhases);
    void setCombination2(const RecurGeom<DIM>& mi1, const RecurGeom<DIM>& mi2,
        map<PhaseType, PhaseType> map_);
    void setMask(const RecurGeom<DIM>& mi1, const RecurGeom<DIM>& mi2, const RecurGeom<DIM>& mask_);
private:
    //! tracks all the phases
    vector<PhaseType> allPhases;
};

}  // namespace merope

#include "RecurMesoStructure.ixx"