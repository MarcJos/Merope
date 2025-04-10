//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"
#include "../../../GenericTools/Variant_shared_ptr.hxx"
#include "../../../Geometry/include/GeomTools.hxx"


namespace merope {

namespace auxiMicroStructure {
//! facilitator to build functions
std::function<PhaseType(PhaseType, PhaseType)> replace_by(
    map<PhaseType, PhaseType>);

enum class TypeOfCombination {
    Void,
    Simple,         // the microstructure is just a simple BasicStruct
    Mask,           // the microstructure is a combination of 2 RecurGeoms, with a function defining the resulting RecurGeom
    Combination2    // the microstructure is a combination of 3 RecurGeoms, the 2nd being a mask
};
}  // namespace auxiMicroStructure

template<unsigned short DIM, class BasicStruct, class BasicType>
class RecursiveStructure {
public:
    using SELF_TYPE = RecursiveStructure<DIM, BasicStruct, BasicType>;
    using BASIC_STRUCT = BasicStruct;
    using VARIANT_TYPE = variant_shared_ptr<BasicStruct,
        tuple<SELF_TYPE, SELF_TYPE, SELF_TYPE>,
        tuple<SELF_TYPE, SELF_TYPE, std::function<BasicType(BasicType, BasicType)>>,
        tuple<SELF_TYPE, SELF_TYPE>
    >;

public:
    //! Default initialization
    RecursiveStructure();
    //! destructor
    virtual ~RecursiveStructure() {}
    //! copy and move constructors. The deep copy is used here.
    //! getter
    const SELF_TYPE& getMask() const;
    template<int Structure_Index>
    const SELF_TYPE& getRecurStructure() const;
    const BasicStruct& getBasicStructure() const;
    //! \return one basic structure which enters the RecursiveStructure
    const BasicStruct& getOneBasicStructure() const;
    std::function<BasicType(BasicType, BasicType)> getTransformFunction() const;
    auxiMicroStructure::TypeOfCombination getTypeOf() const { return typeOfCombination; }
    //! tests whether the structure is not empty and is consistent
    bool isValid() const;

    //protected:
    //! setter
    void setBasicStruct(const BasicStruct& mi);
    void setCombination2(const SELF_TYPE& mi1,
        const SELF_TYPE& mi2,
        std::function<BasicType(BasicType, BasicType)> func);
    void setCombination2(const SELF_TYPE& mi1,
        const SELF_TYPE& mi2,
        map<BasicType, BasicType> map_);
    void setMask(const SELF_TYPE& mi1,
        const SELF_TYPE& mi2,
        const SELF_TYPE& mask_);


private:
    //! type of combination
    auxiMicroStructure::TypeOfCombination typeOfCombination;
    // memory
    VARIANT_TYPE memory_variant;
};

}  // namespace merope


#include "RecurStructure.ixx"


