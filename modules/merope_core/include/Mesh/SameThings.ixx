//! Copyright : see license.txt
//!
//! \brief
//
#ifndef MESH_SAMETHINGS_IXX_
#define MESH_SAMETHINGS_IXX_

#include "../MeropeNamespace.hxx"


namespace merope {
namespace mesh {
namespace sameThings {

template<class DICT_THING, class COMPARISON_FUNCTION>
inline vector<SameThings<Identifier>> getReplacementList(
    const DICT_THING& dictThing, COMPARISON_FUNCTION comparisonFunction) {
    auto vecSameThings = auxi::findSamePairs(dictThing, comparisonFunction);
    return auxi::replaceGraph<typename DICT_THING::mapped_type>(vecSameThings);
}


namespace auxi {
template<class DICT_THING, class COMPARISON_FUNCTION>
inline vector<SameThings<typename DICT_THING::mapped_type>> findSamePairs(const DICT_THING& dictThing, COMPARISON_FUNCTION comparisonFunction) {
    vector<SameThings<typename DICT_THING::mapped_type>> vecSameThings{};
    for (auto it1 = dictThing.begin(); it1 != dictThing.end(); it1++) {
        for (auto it2 = it1; it2 != dictThing.end(); it2++) {
            const auto& x1 = (*it1).second, x2 = (*it2).second;
            AreSame areSame = comparisonFunction(x1, x2);
            if (abs(areSame) == 1 and x1.identifier != x2.identifier) {
                vecSameThings.push_back(make_tuple(x1, x2, areSame));
            }
        }
    }
    return vecSameThings;
}

template<class THING>
inline vector<SameThings<Identifier>> replaceGraph(const vector<SameThings<THING>>& vecSameThings) {
    // Initialization
    //! the map between the idSort to the desired object
    map<IdentifierSort, const THING*> idSort_toObject{};
    //! the graph for going down
    vector<tuple<IdentifierSort, IdentifierSort, AreSame>> idSortObjects{};
    for (size_t i = 0; i < vecSameThings.size(); i++) {
        idSort_toObject[get<0>(vecSameThings[i]).getId_forSort()] = &(get<0>(vecSameThings[i]));
        idSort_toObject[get<1>(vecSameThings[i]).getId_forSort()] = &(get<1>(vecSameThings[i]));
        idSortObjects.push_back(make_tuple(get<0>(vecSameThings[i]).getId_forSort(), get<1>(vecSameThings[i]).getId_forSort(), get<2>(vecSameThings[i])));
    }
    // compute the replacement list
    DictPointer dictPointer(idSortObjects);
    auto vecSameThing_orderedIndex = dictPointer.getOrdered();
    // retrieve the correct replacement list
    vector<SameThings<Identifier>> replacementList{};
    for (const auto& orderedIndex : vecSameThing_orderedIndex) {
#ifndef NDEBUG
        if (idSort_toObject.find(get<0>(orderedIndex)) == idSort_toObject.end() or idSort_toObject.find(get<1>(orderedIndex)) == idSort_toObject.end()) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Problem");
        }
#endif /* NDEBUG */
        replacementList.push_back(make_tuple(idSort_toObject[get<0>(orderedIndex)]->identifier,
            idSort_toObject[get<1>(orderedIndex)]->identifier,
            get<2>(orderedIndex)));
    }
    return replacementList;
}

}// namespace auxi
}// namespace sameThings
}// namespace mesh
} // namespace merope

#endif /* MESH_SAMETHINGS_IXX_ */
