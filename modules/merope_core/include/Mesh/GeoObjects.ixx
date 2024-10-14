//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../MeropeNamespace.hxx"


namespace merope {
namespace mesh {

template<geoObjects::TypeRelation RELATION>
inline geoObjects::GeoObject::TYPE<RELATION>& geoObjects::GeoObject::get() {
    if constexpr (RELATION == TypeRelation::Root) { return roots; };
    if constexpr (RELATION == TypeRelation::Leaf) { return leaves; };
}

template<geoObjects::TypeRelation RELATION>
inline const geoObjects::GeoObject::TYPE<RELATION>& geoObjects::GeoObject::get() const {
    if constexpr (RELATION == TypeRelation::Root) { return roots; };
    if constexpr (RELATION == TypeRelation::Leaf) { return leaves; };
}

template<geoObjects::TypeRelation RELATION>
inline vector<Identifier> geoObjects::GeoObject::find(Identifier id) const {
    auto& myList = get<RELATION>();
    vector<Identifier> result{};
    for (auto it = myList.begin(); it != myList.end(); it++) {
        if (abs(*it) == id) {
            result.push_back(*it);
        }
    }
    return result;
}

//! PerSurface
template<unsigned short DIM>
inline void geoObjects::PerSurface<DIM>::swapSurf() {
    if (not isSingular()) {
        leaves = { leaves[1], leaves[0] };
        translation = -translation;
    }
}

//!
template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES>
inline bool geoObjects::findSingularAndRemove(DICT_ROOTS& dictRoots, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves) {
    vector<Identifier> singularThings = findSingularThings(dictThings);
    removeSingularThings(singularThings, dictRoots, dictThings, dictLeaves);
    return singularThings.size() > 0;
}



template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES,
    class DICT_PER_THINGS>
inline bool mesh::geoObjects::findSingularAndRemovePeriodic(DICT_ROOTS& dictRoots, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves, DICT_PER_THINGS& dictPerThings) {
    vector<Identifier> singularThings = findSingularThings(dictThings);
    removePeriodicRoots(singularThings, dictPerThings, dictThings);
    removeSingularThings(singularThings, dictRoots, dictThings, dictLeaves);
    return singularThings.size() > 0;
}

template<class DICT_THINGS>
inline vector<Identifier> geoObjects::findSingularThings(const DICT_THINGS& dictThings) {
    vector<Identifier> singularThings = {};
    for (auto it = dictThings.begin(); it != dictThings.end(); it++) {
        if (it->second.isSingular()) {
            singularThings.push_back(it->first);
        }
    }
    return singularThings;
}

template<class DICT_PER_THINGS, class DICT_THINGS>
inline void mesh::geoObjects::removePeriodicRoots(const vector<Identifier>& singularThings, DICT_PER_THINGS& dictPerThings, const DICT_THINGS& dictThings) {
    for (auto id : singularThings) {
        const auto& thing = dictThings.at(id);
        if (thing.isPeriodic()) {
            auto id_periodic_root = thing.getPeriodicRoot();
            dictPerThings.at(id_periodic_root).removeLeafSingular(id);
            if (dictPerThings.at(id_periodic_root).isSingular()) {
                dictPerThings.erase(id_periodic_root);
            }
        }
    }
}


template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES>
inline void geoObjects::removeSingularThings(const vector<Identifier>& singularThings, DICT_ROOTS& dictRoots,
    DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves) {
    for (auto identifier : singularThings) {
        auto thing = dictThings.at(identifier);
        for (const auto& root_id : thing.template get<TypeRelation::Root>()) {
            dictRoots.at(root_id).removeLeafSingular(thing.identifier);
        }
        removeObjectDownwards(thing.identifier, dictThings, dictLeaves);
    }
}

template<class DICT_THINGS, class DICT_LEAVES>
inline void geoObjects::removeObjectDownwards(Identifier idObject, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves) {
    for (const auto& id_leave : dictThings.at(idObject).template get<TypeRelation::Leaf>()) {
        auto& leaf = dictLeaves.at(abs(id_leave));
        leaf.removeRoot(idObject);
    }
    dictThings.erase(idObject);  // suppress the thing to be replaced
}


template<class OBJ>
inline geoObjects::AreSame geoObjects::areSame(const OBJ& obj1, const OBJ& obj2) {
    //! Edge or PerSurface
    if constexpr (std::is_same_v<OBJ, Edge> or std::is_same_v<OBJ, PerSurface<3>>) {
        // DEBUG
        obj1.bugTest();
        obj2.bugTest();
        // DEBUG
        Identifier x0 = obj1.leaves[0], x1 = obj1.leaves[1], y0 = obj2.leaves[0], y1 = obj2.leaves[1];
        if (x0 == y0 and x1 == y1) {
            return AreSame::Same;
        } else if (x0 == y1 and x1 == y0) {
            return AreSame::Reverse;
        } else {
            return AreSame::Different;
        }
    }
    //! CurveLoop, Surface, SurfaceLoop, Solid
    else if constexpr (std::is_same_v<OBJ, CurveLoop>
        or std::is_same_v<OBJ, Surface>
        or std::is_same_v<OBJ, SurfaceLoop>
        or std::is_same_v<OBJ, Solid>) {
        auto leaves1 = obj1.leaves, leaves2 = obj2.leaves;
        std::sort(leaves1.begin(), leaves1.end());
        std::sort(leaves2.begin(), leaves2.end());
        if (leaves1 == leaves2) {
            return AreSame::Same;
        }
        for (auto& leaf : leaves2) {
            leaf *= -1;
        }
        std::sort(leaves2.begin(), leaves2.end());
        if (leaves1 == leaves2) {
            return AreSame::Reverse;
        }
        return AreSame::Different;
    }
    //! PerPoint
    else if constexpr (std::is_same_v<OBJ, PerPoint>) {
        return static_cast<AreSame>(std::find_first_of(obj1.leaves.begin(), obj1.leaves.end(), obj2.leaves.begin(), obj2.leaves.end()) != obj1.leaves.end());
    }
    return AreSame::Different;
}


//!
template<class DICT_THING, class DICT_LEAVES>
inline bool geoObjects::isGraphCoherent(const DICT_THING& dict_thing,
    const DICT_LEAVES& dict_leaves) {
    return isGraphCoherent_auxi<TypeRelation::Leaf>(dict_thing, dict_leaves) and isGraphCoherent_auxi<TypeRelation::Root>(dict_leaves, dict_thing);
}

template<geoObjects::TypeRelation RELAT, class DICT_THING, class DICT_RELATED>
inline bool geoObjects::isGraphCoherent_auxi(const DICT_THING& dict_thing, const DICT_RELATED& dict_related) {
    for (const auto& elem : dict_thing) {
        const auto& id_thing = elem.first;
        const auto& thing = elem.second;
        for (auto index_thing_related : thing.template get<RELAT>()) {
            index_thing_related = abs(index_thing_related);
            auto found = dict_related.find(index_thing_related);
            if (found == dict_related.end() or found->second.template find<Invert(RELAT)>(id_thing).size() == 0) {
                cerr << getName(thing.name) << " " << id_thing << endl;
                return false;
            }
        }
    }
    return true;
}

//! merge
template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES>
inline void geoObjects::vec_merge(vector<SameThings<Identifier>> vecSameThings, DICT_ROOTS& dictRoots, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves) {
    for (const auto& sameThings : vecSameThings) {
        merge(sameThings, dictRoots, dictThings, dictLeaves);
    }
}


template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES>
void geoObjects::merge(const SameThings<Identifier> sameThings, DICT_ROOTS& dictRoots, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves) {
    auto thing_original = dictThings.at(get<0>(sameThings));
    auto& thing_replacer = dictThings.at(get<1>(sameThings));
    auto sense = get<2>(sameThings);
    //! replace the thing for the roots
    for (const auto& id_root : thing_original.template get<TypeRelation::Root>()) {
        auto& root = dictRoots.at(id_root);
        root.replaceLeaf(thing_original.identifier, thing_replacer.identifier, sense);
        thing_replacer.addRoot(id_root);
    }
    //! replace the thing for the leaves
    for (const auto& id_leaf : thing_original.template get<TypeRelation::Leaf>()) {
        auto& leaf = dictLeaves.at(abs(id_leaf));
        leaf.removeRoot(thing_original.identifier);
        leaf.addRoot(thing_replacer.identifier);
    }
    //! removes the object
    dictThings.erase(thing_original.identifier);
}

template<class DICT_THING>
inline vector<sameThings::SameThings<Identifier>> sameThings::getReplacementList(const DICT_THING& dictThing) {
    auto comparisonFunction = [](const auto& obj1, const auto& obj2) {return geoObjects::areSame(obj1, obj2);};
    return getReplacementList(dictThing, comparisonFunction);
}

template<class DICT_THINGS>
inline vector<typename DICT_THINGS::mapped_type> geoObjects::findSingular(
    DICT_THINGS dictThings) {
    vector<typename DICT_THINGS::mapped_type> vecSingularThings{};
    for (const auto& thing : dictThings) {
        if (thing.second.isSingular()) {
            vecSingularThings.push_back(thing.second);
        }
    }
    return vecSingularThings;
}

}  // namespace mesh
}  // namespace merope



