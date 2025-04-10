//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace merope {

template<unsigned short DIM, class BasicStruct, class BasicType>
RecursiveStructure<DIM, BasicStruct, BasicType>::RecursiveStructure() :
    typeOfCombination{ auxiMicroStructure::TypeOfCombination::Void },
    memory_variant{ VARIANT_TYPE(shared_ptr<BasicStruct>(nullptr)) } {}

template<unsigned short DIM, class BasicStruct, class BasicType>
void RecursiveStructure<DIM, BasicStruct, BasicType>::setBasicStruct(
    const BasicStruct& basicStruct_) {
    typeOfCombination = auxiMicroStructure::TypeOfCombination::Simple;
    memory_variant = shared_ptr<BasicStruct>(new BasicStruct(basicStruct_));
    isValid();
}

template<unsigned short DIM, class BasicStruct, class BasicType>
void RecursiveStructure<DIM, BasicStruct, BasicType>::setCombination2(
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct1,
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct2,
    std::function<BasicType(BasicType, BasicType)> func) {
    typeOfCombination = auxiMicroStructure::TypeOfCombination::Combination2;
    memory_variant = shared_ptr<tuple<SELF_TYPE, SELF_TYPE, std::function<BasicType(BasicType, BasicType)>>>(
        new tuple<SELF_TYPE, SELF_TYPE, std::function<BasicType(BasicType, BasicType)>>{
            basicStruct1, basicStruct2, func });
    isValid();
}

template<unsigned short DIM, class BasicStruct, class BasicType>
void RecursiveStructure<DIM, BasicStruct, BasicType>::setCombination2(
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct1,
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct2,
    map<BasicType, BasicType> map_) {
    setCombination2(basicStruct1, basicStruct2, auxiMicroStructure::replace_by(map_));
    isValid();
}

template<unsigned short DIM, class BasicStruct, class BasicType>
void RecursiveStructure<DIM, BasicStruct, BasicType>::setMask(
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct1,
    const RecursiveStructure<DIM, BasicStruct, BasicType>& basicStruct2,
    const RecursiveStructure<DIM, BasicStruct, BasicType>& mask_) {
    typeOfCombination = auxiMicroStructure::TypeOfCombination::Mask;
    memory_variant = shared_ptr<tuple<SELF_TYPE, SELF_TYPE, SELF_TYPE>>(
        new tuple<SELF_TYPE, SELF_TYPE, SELF_TYPE>{ basicStruct1, basicStruct2, mask_ });
    isValid();
}

template<unsigned short DIM, class BasicStruct, class BasicType>
const BasicStruct& RecursiveStructure<DIM, BasicStruct, BasicType>::getBasicStructure() const {
    return *(std::get<0>(memory_variant));
}

template<unsigned short DIM, class BasicStruct, class BasicType>
const BasicStruct& RecursiveStructure<DIM, BasicStruct, BasicType>::getOneBasicStructure() const {
    switch (this->getTypeOf()) {
    case auxiMicroStructure::TypeOfCombination::Simple:
        return getBasicStructure();
    default:
        return getRecurStructure<0>().getOneBasicStructure();
    }
}

template<unsigned short DIM, class BasicStruct, class BasicType>
const RecursiveStructure<DIM, BasicStruct, BasicType>&
RecursiveStructure<DIM, BasicStruct, BasicType>::getMask() const {
    Merope_assert(this->getTypeOf() == auxiMicroStructure::TypeOfCombination::Mask, "Incorrect type of Combination");
    return std::get<2>(*(std::get<1>(memory_variant)));
}

template<unsigned short DIM, class BasicStruct, class BasicType>
template<int Structure_Index>
const RecursiveStructure<DIM, BasicStruct, BasicType>& RecursiveStructure<DIM, BasicStruct, BasicType>::getRecurStructure() const {
    switch (this->getTypeOf()) {
    case auxiMicroStructure::TypeOfCombination::Mask:
        return std::get<Structure_Index>(*(std::get<1>(memory_variant)));
    case auxiMicroStructure::TypeOfCombination::Combination2:
        return std::get<Structure_Index>(*(std::get<2>(memory_variant)));
    default:
        Merope_assert(false, "Incorrect type of Combination");
    }
}

template<unsigned short DIM, class BasicStruct, class BasicType>
std::function<BasicType(BasicType, BasicType)> RecursiveStructure<DIM,
    BasicStruct, BasicType>::getTransformFunction() const {
    Merope_assert(this->getTypeOf() == auxiMicroStructure::TypeOfCombination::Combination2, "Incorrect type of Combination");
    return std::get<2>(*(std::get<2>(memory_variant)));
}

template<unsigned short DIM, class BasicStruct, class BasicType>
bool RecursiveStructure<DIM, BasicStruct, BasicType>::isValid() const {
    switch (this->getTypeOf()) {
    case auxiMicroStructure::TypeOfCombination::Simple:
    {
        return std::get_if<shared_ptr<BasicStruct>>(&memory_variant);
    }
    case auxiMicroStructure::TypeOfCombination::Combination2:
    {
        return areCompatible(getRecurStructure<0>(), getRecurStructure<1>());
    }
    case auxiMicroStructure::TypeOfCombination::Mask:
    {
        return areCompatible(getRecurStructure<0>(), getRecurStructure<1>()) and areCompatible(getRecurStructure<0>(), getMask());
    }
    default:
    {
        Merope_assert(false, "Not valid");
    }
    }
    return true;
}

}  // namespace merope


