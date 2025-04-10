//! Copyright : see license.txt
//!
//! \brief

#pragma once

namespace merope {
template<unsigned short DIM>
RecurGeom<DIM>::RecurGeom() :
    RecursiveStructure_inTorus<DIM, BasicStruct, BasicType>(),
    allPhases{}{
}

template<unsigned short DIM>
void RecurGeom<DIM>::setMultiInclusions(const MultiInclusions<DIM>& mi) {
    RecursiveStructure<DIM, BasicStruct, BasicType>::setBasicStruct(mi);
    allPhases = mi.getAllPhases();
}

template<unsigned short DIM>
void RecurGeom<DIM>::setCombination2(const RecurGeom<DIM>& mi1,
    const RecurGeom<DIM>& mi2,
    std::function<PhaseType(PhaseType, PhaseType)> func,
    vector<PhaseType> allPhases_) {
    RecursiveStructure<DIM, BasicStruct, BasicType>::setCombination2(mi1, mi2, func);
    allPhases = allPhases_;
}

template<unsigned short DIM>
void RecurGeom<DIM>::setCombination2(const RecurGeom<DIM>& mi1,
    const RecurGeom<DIM>& mi2, map<PhaseType, PhaseType> map_) {
    RecursiveStructure<DIM, BasicStruct, BasicType>::setCombination2(mi1, mi2, map_);
    auto allPhases_ = mi1.getAllPhases();
    for (const auto& [key, value] : map_) {
        allPhases_.push_back(value);
    }
    allPhases = allPhases_;
}

template<unsigned short DIM>
void RecurGeom<DIM>::setMask(const RecurGeom<DIM>& mi1,
    const RecurGeom<DIM>& mi2, const RecurGeom<DIM>& mask_) {
    RecursiveStructure<DIM, BasicStruct, BasicType>::setMask(mi1, mi2, mask_);
    allPhases = mi1.getAllPhases();
    for (const auto& ph : mi2.getAllPhases()) {
        allPhases.push_back(ph);
    }
}
}  // namespace  merope
