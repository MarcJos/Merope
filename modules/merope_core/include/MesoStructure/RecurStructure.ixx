//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef PREMICROSTRUCTURE_IXX_
#define PREMICROSTRUCTURE_IXX_

#include "../Grid/AreCompatible.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM, class BasicStruct, class BasicType>
inline RecursiveStructure<DIM, BasicStruct, BasicType>::RecursiveStructure(
    auxiMicroStructure::TypeOfCombination typeOf_):
    InsideTorus<DIM>(),
    typeOfCombination{ typeOf_ },
    basicStructure{ nullptr },
    recurStructure1{ nullptr },
    recurStructure2{ nullptr },
    transformFunction{ nullptr },
    mask{ nullptr }{
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline RecursiveStructure<DIM, BasicStruct, BasicType>::RecursiveStructure(
    const RecursiveStructure<DIM, BasicStruct, BasicType>& other)
    : RecursiveStructure<DIM, BasicStruct, BasicType>(other.getTypeOf()) {
    switch (other.getTypeOf()) {
    case auxiMicroStructure::TypeOfCombination::Void:
    {
        break;
    }
    case auxiMicroStructure::TypeOfCombination::Simple:
    {
        setBasicStruct(other.getBasicStructure());
        break;
    }
    case auxiMicroStructure::TypeOfCombination::Combination2:
    {
        setCombination2(other.getRecurStructure1(), other.getRecurStructure2(), other.getTransformFunction());
        break;
    }
    case auxiMicroStructure::TypeOfCombination::Mask:
    {
        setMask(other.getRecurStructure1(), other.getRecurStructure2(), other.getMask());
        break;
    }
    default:
    {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Not programmed");
    }
    }
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline RecursiveStructure<DIM, BasicStruct, BasicType>& RecursiveStructure<DIM,
    BasicStruct, BasicType>::operator =(
        const RecursiveStructure<DIM, BasicStruct, BasicType>& other) {
    RecursiveStructure<DIM, BasicStruct, BasicType> that(other);
    swap(that);
    return *this;
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline RecursiveStructure<DIM, BasicStruct, BasicType>& RecursiveStructure<DIM,
    BasicStruct, BasicType>::operator =(
        RecursiveStructure<DIM, BasicStruct, BasicType>&& other) {
    RecursiveStructure<DIM, BasicStruct, BasicType> that(std::move(other));
    swap(that);
    return *this;
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline void RecursiveStructure<DIM, BasicStruct, BasicType>::swap(
    RecursiveStructure<DIM, BasicStruct, BasicType>& other) {
    std::swap(this->typeOfCombination, other.typeOfCombination);
    std::swap(this->basicStructure, other.basicStructure);
    std::swap(this->recurStructure1, other.recurStructure1);
    std::swap(this->recurStructure2, other.recurStructure2);
    std::swap(this->transformFunction, other.transformFunction);
    std::swap(this->mask, other.mask);
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline void RecursiveStructure<DIM, BasicStruct, BasicType>::setBasicStruct(
    const BasicStruct& basicStruct_) {
    this->setLength(basicStruct_.getL());
    typeOfCombination = auxiMicroStructure::TypeOfCombination::Simple;
    basicStructure.reset(new BasicStruct(basicStruct_));
    isValid();
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline void RecursiveStructure<DIM, BasicStruct, BasicType>::setCombination2(
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct1,
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct2,
    std::function<BasicType(BasicType, BasicType)> func) {
    this->setLength(basicStruct1.getL());
    typeOfCombination = auxiMicroStructure::TypeOfCombination::Combination2;
    recurStructure1.reset(new RecursiveStructure<DIM, BasicStruct, BasicType>(basicStruct1));
    recurStructure2.reset(new RecursiveStructure<DIM, BasicStruct, BasicType>(basicStruct2));
    transformFunction.reset(new std::function<BasicType(BasicType, BasicType)>(func));
    isValid();
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline void RecursiveStructure<DIM, BasicStruct, BasicType>::setCombination2(
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct1,
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct2,
    map<BasicType, BasicType> map_) {
    setCombination2(basicStruct1, basicStruct2, auxiMicroStructure::replace_by(map_));
    isValid();
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline void RecursiveStructure<DIM, BasicStruct, BasicType>::setMask(
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct1,
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct2,
    const RecursiveStructure<DIM, BasicStruct, BasicType>& mask_) {
    this->setLength(basicStruct1.getL());
    typeOfCombination = auxiMicroStructure::TypeOfCombination::Mask;
    recurStructure1.reset(new RecursiveStructure<DIM, BasicStruct, BasicType>(basicStruct1));
    recurStructure2.reset(new RecursiveStructure<DIM, BasicStruct, BasicType>(basicStruct2));
    mask.reset(new RecursiveStructure<DIM, BasicStruct, BasicType>(mask_));
    isValid();
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline const BasicStruct& RecursiveStructure<DIM, BasicStruct, BasicType>::getBasicStructure() const {
    if (basicStructure) {
        return *basicStructure;
    }
    else {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline const BasicStruct& RecursiveStructure<DIM, BasicStruct, BasicType>::getOneBasicStructure() const {
    if (basicStructure) {
        return *basicStructure;
    }
    else {
        return getRecurStructure1().getOneBasicStructure();
    }
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline const RecursiveStructure<DIM, BasicStruct, BasicType>&
RecursiveStructure<DIM, BasicStruct, BasicType>::getMask() const {
    if (mask) {
        return *mask;
    }
    else {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline const RecursiveStructure<DIM, BasicStruct, BasicType>&
RecursiveStructure<DIM, BasicStruct, BasicType>::getRecurStructure1() const {
    if (recurStructure1) {
        return *recurStructure1;
    }
    else {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline const RecursiveStructure<DIM, BasicStruct, BasicType>&
RecursiveStructure<DIM, BasicStruct, BasicType>::getRecurStructure2() const {
    if (recurStructure2) {
        return *recurStructure2;
    }
    else {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline std::function<BasicType(BasicType, BasicType)> RecursiveStructure<DIM,
    BasicStruct, BasicType>::getTransformFunction() const {
    if (transformFunction) {
        return *transformFunction;
    }
    else {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM, class BasicStruct, class BasicType>
inline bool RecursiveStructure<DIM, BasicStruct, BasicType>::isValid() const {
    switch (this->getTypeOf()) {
    case auxiMicroStructure::TypeOfCombination::Simple:
    {
        return true;
    }
    case auxiMicroStructure::TypeOfCombination::Combination2:
    {
        return areCompatible(getRecurStructure1(), getRecurStructure2());
    }
    case auxiMicroStructure::TypeOfCombination::Mask:
    {
        return areCompatible(getRecurStructure1(), getRecurStructure2()) and areCompatible(getRecurStructure1(), getMask());
    }
    default:
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

////

template<unsigned short DIM>
inline RecurGeom<DIM>::RecurGeom():
    RecursiveStructure<DIM, BasicStruct, BasicType>(),
    allPhases{}{
}

template<unsigned short DIM>
inline void RecurGeom<DIM>::setMultiInclusions(const MultiInclusions<DIM>& mi) {
    RecursiveStructure<DIM, BasicStruct, BasicType>::setBasicStruct(mi);
    allPhases = mi.getAllPhases();
}

template<unsigned short DIM>
inline void RecurGeom<DIM>::setCombination2(const RecurGeom<DIM>& mi1,
    const RecurGeom<DIM>& mi2,
    std::function<PhaseType(PhaseType, PhaseType)> func,
    vector<PhaseType> allPhases_) {
    RecursiveStructure<DIM, BasicStruct, BasicType>::setCombination2(mi1, mi2, func);
    allPhases = allPhases_;
}

template<unsigned short DIM>
inline void RecurGeom<DIM>::setCombination2(const RecurGeom<DIM>& mi1,
    const RecurGeom<DIM>& mi2, map<PhaseType, PhaseType> map_) {
    RecursiveStructure<DIM, BasicStruct, BasicType>::setCombination2(mi1, mi2, map_);
    auto allPhases_ = mi1.getAllPhases();
    for (const auto& [key, value] : map_) {
        allPhases_.push_back(value);
    }
    allPhases = allPhases_;
}

template<unsigned short DIM>
inline void RecurGeom<DIM>::setMask(const RecurGeom<DIM>& mi1,
    const RecurGeom<DIM>& mi2, const RecurGeom<DIM>& mask_) {
    RecursiveStructure<DIM, BasicStruct, BasicType>::setMask(mi1, mi2, mask_);
    allPhases = mi1.getAllPhases();
    for (const auto& ph : mi2.getAllPhases()) {
        allPhases.push_back(ph);
    }
}

} // namespace merope

#endif /* PREMICROSTRUCTURE_IXX_ */
