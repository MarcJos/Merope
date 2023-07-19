//! Copyright : see license.txt
//!
//! \brief
//
#include "Mesh/SameThings.hxx"
#include "../../AlgoPacking/src/AuxiFunctions.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

namespace mesh {
namespace sameThings {
namespace auxi {

//! DictPointer

void DictPointer::Initialize() {
    //! initialize
    for (const auto& sameThings : vecSameThings) {
        const auto& x0 = get<0>(sameThings), x1 = get<1>(sameThings);
        internalGraph[x0] = {};
        internalGraph[x1] = {};
    }
    //! fill
    for (const auto& sameThings : vecSameThings) {
        const auto& x0 = get<0>(sameThings), x1 = get<1>(sameThings);
        auto areSame = get<2>(sameThings);
        internalGraph[x0].push_back(make_tuple(x1, areSame));
        internalGraph[x1].push_back(make_tuple(x0, areSame));
    }
}

bool DictPointer::EnrichGraph() {
    bool hasChanged = false;
    for (auto it = internalGraph.begin(); it != internalGraph.end(); it++) {
        auto& relatedElements = it->second;
        for (size_t i = 0; i < relatedElements.size(); i++) {
            const auto& elem = relatedElements[i];
            const auto& relatedElements2 = internalGraph[get<0>(elem)];
            for (const auto& elem2 : relatedElements2) {
                bool elem2_notContainedIn_relatedElements = find_if(relatedElements.begin(), relatedElements.end(), [elem2](const auto& e) {
                    return get<0>(e) == get<0>(elem2);
                    }) == relatedElements.end();
                    if (elem2_notContainedIn_relatedElements) {
                        auto newId = get<0>(elem2);
                        auto newAreSame = static_cast<AreSame>(get<1>(elem) * get<1>(elem2));
                        relatedElements.push_back(make_tuple(newId, newAreSame));
                        hasChanged = true;
                    }
            }
        }
    }
    return hasChanged;
}

vector<mesh::sameThings::SameThings<IdentifierSort> > DictPointer::getOrdered() {
    this->Initialize();
    while (this->EnrichGraph()) {} // fixme : should work if done only once
    return getResult();
}

vector<SameThings<IdentifierSort>> DictPointer::getResult() const {
    vector<SameThings<IdentifierSort>> vecSameThing{};
    for (auto it = internalGraph.begin(); it != internalGraph.end(); it++) {
        auto id = it->first;
        const auto& vecRelated = it->second;
        auto replacement = *min_element(vecRelated.begin(), vecRelated.end(), [](const auto& duet1, const auto& duet2) {
            return get<0>(duet1) < get<0>(duet2);
            });
        vecSameThing.push_back(make_tuple(id, get<0>(replacement), get<1>(replacement)));
    }
    auxi_function::erase_if(vecSameThing, [](const auto& triplet) {return get<0>(triplet) == get<1>(triplet);});
    return vecSameThing;
}

} // namespace auxi
} // namespace sameThings
} // namespace mesh
} // namespace merope



