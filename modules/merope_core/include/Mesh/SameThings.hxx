//! Copyright : see license.txt
//!
//! \briefFor merging equal or almost-equal geometrical objects
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/Geometry/GeomTypes.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

namespace mesh {
namespace sameThings {

//! to capture the relationship between two objects
enum AreSame {
    Same = 1,           // same object, same orientation
    Reverse = -1,       // same object, reverse orientation
    Different = 0       // different object
};

//! to capture the relationship between two objects
template<class C>
using SameThings = tuple<C, C, AreSame>;

//! define an order on the geometrical objects.
//! the object with the larger IdentifierSort should be merged into the object with the lesser one
struct IdentifierSort : public array<Identifier, 2>{
    IdentifierSort(const array<Identifier, 2>& array_) : array<Identifier, 2>(array_) {}
    bool operator<(IdentifierSort ids2) { return std::lexicographical_compare(this->begin(), this->end(), ids2.begin(), ids2.end()); }
    bool operator>(IdentifierSort ids2) { return std::lexicographical_compare(this->begin(), this->end(), ids2.begin(), ids2.end(), [](const auto& o1, const auto& o2) {return o1 > o2;}); }
    bool operator==(IdentifierSort ids2) { return (*this)[0] == ids2[0] and (*this)[1] == ids2[1]; }
};

inline ostream& operator<<(ostream& os, const IdentifierSort& is) {
    os << "(" << is[0] << ',' << is[1] << ')';
    return os;
}

//! \return same entities, according to a comparison function
//! \param dictThing : map of things
//! \param comparisonFunction : function to compare two entities
template<class DICT_THING, class COMPARISON_FUNCTION>
vector<SameThings<Identifier>> getReplacementList(const DICT_THING& dictThing, COMPARISON_FUNCTION comparisonFunction);

/////////////////////////
//! Auxiliary functions
/////////////////////////
namespace auxi {

//! \return same entities, according to a comparison function
//! \param dictThing : map of things
//! \param comparisonFunction : function to compare two entities
template<class DICT_THING, class COMPARISON_FUNCTION>
vector<SameThings<typename DICT_THING::mapped_type>> findSamePairs(const DICT_THING& dictThing, COMPARISON_FUNCTION comparisonFunction);

//! \return a list of replacements, from a list of same entities to , where all the same entities are to be merged into a single final one
//! \param vecSameThings : vector of things that are the same
template<class THING>
vector<SameThings<Identifier>> replaceGraph(const vector<SameThings<THING>>& vecSameThings);

//! class for turning a vector of same objects into a vector of replacements
struct DictPointer {
public:
    //! constructor
    DictPointer(const vector<mesh::sameThings::SameThings<IdentifierSort>>& vecSameThings_) : vecSameThings(vecSameThings_), internalGraph{} {};
    //! apply
    vector<mesh::sameThings::SameThings<IdentifierSort>> getOrdered();
private:
    //! internal copy of vecSameThings
    const vector<mesh::sameThings::SameThings<IdentifierSort>> vecSameThings;
    //! internal graph
    std::map<IdentifierSort, vector<tuple<IdentifierSort, AreSame>>> internalGraph;
    //! Initialize the graph
    void Initialize();
    //! Enriches the graph
    bool EnrichGraph();
    //! \return the desired replacement list
    vector<SameThings<IdentifierSort>> getResult() const;
};

}  // namespace auxi
}  // namespace sameThings
}  // namespace mesh
}  // namespace merope

#include "../Mesh/SameThings.ixx"

