//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef PREMICROSTRUCTURE_HXX_
#define PREMICROSTRUCTURE_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Field/CartesianField.hxx"
#include "../Geometry/GeomTools.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"


#include "../MeropeNamespace.hxx"


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
} // namespace auxiMicroStructure


template<unsigned short DIM, class BasicStruct, class BasicType>
class RecursiveStructure: public InsideTorus<DIM> {
    //
    static_assert(std::is_base_of<InsideTorus<DIM>, BasicStruct>::value);
    //
public:
    //! Default initialization
    RecursiveStructure(auxiMicroStructure::TypeOfCombination typeOf_ = auxiMicroStructure::TypeOfCombination::Void);
    //! destructor
    virtual ~RecursiveStructure() {};
    //! copy and move constructors. The deep copy is used here.
    RecursiveStructure(const RecursiveStructure<DIM, BasicStruct, BasicType>& other);
    RecursiveStructure(RecursiveStructure<DIM, BasicStruct, BasicType>&& other) = default;
    RecursiveStructure<DIM, BasicStruct, BasicType>& operator=(const RecursiveStructure<DIM, BasicStruct, BasicType>& other);
    RecursiveStructure<DIM, BasicStruct, BasicType>& operator=(RecursiveStructure<DIM, BasicStruct, BasicType>&& other);
    void swap(RecursiveStructure<DIM, BasicStruct, BasicType>& other);
    //! getter
    const RecursiveStructure<DIM, BasicStruct, BasicType>& getMask() const;
    const RecursiveStructure<DIM, BasicStruct, BasicType>& getRecurStructure1() const;
    const RecursiveStructure<DIM, BasicStruct, BasicType>& getRecurStructure2() const;
    const BasicStruct& getBasicStructure() const;
    //! \return one basic structure which enters the RecursiveStructure
    const BasicStruct& getOneBasicStructure() const;
    std::function<BasicType(BasicType, BasicType)> getTransformFunction() const;
    auxiMicroStructure::TypeOfCombination getTypeOf() const { return typeOfCombination; }
    //! tests whether the structure is not empty and is consistent
    bool isValid() const;

protected:
    //! setter
    void setBasicStruct(const BasicStruct& mi);
    void setCombination2(const RecursiveStructure<DIM, BasicStruct, BasicType>& mi1,
        const RecursiveStructure<DIM, BasicStruct, BasicType>& mi2,
        std::function<BasicType(BasicType, BasicType)> func);
    void setCombination2(const RecursiveStructure<DIM, BasicStruct, BasicType>& mi1,
        const RecursiveStructure<DIM, BasicStruct, BasicType>& mi2,
        map<BasicType, BasicType> map_);
    void setMask(const RecursiveStructure<DIM, BasicStruct, BasicType>& mi1,
        const RecursiveStructure<DIM, BasicStruct, BasicType>& mi2,
        const RecursiveStructure<DIM, BasicStruct, BasicType>& mask_);

private:
    //! type of combination
    auxiMicroStructure::TypeOfCombination typeOfCombination;
    //! case 1 : MultiInclusions
    unique_ptr<BasicStruct> basicStructure;
    //! case 2 : Combination2
    unique_ptr<RecursiveStructure<DIM, BasicStruct, BasicType>> recurStructure1;
    unique_ptr<RecursiveStructure<DIM, BasicStruct, BasicType>> recurStructure2;
    unique_ptr<std::function<BasicType(BasicType, BasicType)>> transformFunction;
    //! case 3 : Mask
    unique_ptr<RecursiveStructure<DIM, BasicStruct, BasicType>> mask;
};

template<unsigned short DIM>
using RecurField = RecursiveStructure<DIM, CartesianField<DIM>, double>;

template<unsigned short DIM>
class RecurGeom: public RecursiveStructure<DIM, MultiInclusions<DIM>, PhaseType> {
    using BasicStruct = MultiInclusions<DIM>;
    using BasicType = PhaseType;
public:
    //! Default initialization
    RecurGeom();
    //! get all contained phases
    const vector<PhaseType>& getAllPhases() const { return allPhases; };
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

} // namespace merope


#include "RecurStructure.ixx"

#endif /* PREMICROSTRUCTURE_HXX_ */
