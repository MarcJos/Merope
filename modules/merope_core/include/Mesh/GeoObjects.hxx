//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"
#include "../Mesh/SameThings.hxx"


namespace merope {
namespace mesh {
namespace geoObjects {
using namespace mesh::sameThings;

//! types of geometrical mesh objects
//! \warning to extend this list, also modify the function mesh::geoObjects::getName
enum class TypeObject {
    Unknown = 0, Point = 1, Edge = 2, CurveLoop = 3, Surface = 4, SurfaceLoop = 5, Solid = 6, PhysicalSurface = 7, PhysicalSolid = 8, PerPoint = 9, PeriodicSurface = 10,
};

//! converts the name of a geometrical mesh object into a string
string getName(TypeObject typeObject);

//! graph relation
enum class TypeRelation {
    Leaf, Root
};

//! invert the graph relation
constexpr inline TypeRelation Invert(TypeRelation RELAT) {
    if (RELAT == TypeRelation::Leaf) return TypeRelation::Root;
    else                            return TypeRelation::Leaf;
}


//! base class for geometrical mesh object
class GeoObject {
    template<TypeRelation RELAT>
    using TYPE = typename std::conditional_t<RELAT == TypeRelation::Leaf, vector<Identifier>, std::set<Identifier>>;
public:
    //! identifier of the object
    Identifier identifier;
    //! related objects of lesser dimension
    TYPE<TypeRelation::Leaf> leaves;
    //! name of the oject
    TypeObject name;
private:
    //! related objects of higher dimension
    TYPE<TypeRelation::Root> roots;
public:
    //! constructor
    GeoObject(TypeObject name_, Identifier identifier_, const vector<Identifier>& leaves_) :
        identifier{ identifier_ }, leaves(leaves_), name{ name_ } {}
    //! destructor
    virtual ~GeoObject() {}
    //! remove all roots
    void resetRoots() { roots = {}; }
    //! remove a root
    void removeRoot(Identifier id) { roots.erase(id); }
    //! removes a leaf
    //! \warning : only if the leaf is singular. Otherwise, unexpected behaviour
    void removeLeafSingular(Identifier id);
    //! add a root
    void addRoot(Identifier id) { roots.insert(id); }
    //! replace a leaf
    void replaceLeaf(Identifier id_1, Identifier id_2, long sense);
    //! return the identifier for sorting (when merging)
    virtual IdentifierSort getId_forSort() const { return IdentifierSort({ 0, identifier }); }
    //! test whether the object is cohrent or not
    virtual bool isCoherent()const { return true; }
    //! throws an error in case of problem
    //! fixme : when it works, remove it
    void bugTest() const { if (not isCoherent()) throw runtime_error("The object " + getName(name) + " is not coherent"); }
    //! test whether the object is singular or not
    virtual bool isSingular()const { bugTest(); return leaves.size() == 0; }
    //! return the desired list
    template<TypeRelation RELAT>
    TYPE<RELAT>& get();
    template<TypeRelation RELAT>
    const TYPE<RELAT>& get() const;
    //! print the thing
    virtual void print(std::ostream& f) const;
    //! shift the indices by a given shift
    virtual void shiftIndices(const Identifier shift);
    //! find a leaf/root
    template<TypeRelation>
    vector<Identifier> find(Identifier id) const;
};

//! base class for geometrical mesh object that connects periodic copies
class PerGeoObject : public GeoObject {
public:
    //! constructor
    PerGeoObject(TypeObject name_, Identifier identifier_, vector<Identifier> leaves_) : GeoObject(name_, identifier_, leaves_) {}
    //! add a leaf
    void addLeaf(Identifier id) { leaves.push_back(id); }
};

//! base class for geometrical mesh object the periodic copies of which are tracked by a PerGeoObject
class GeoObject_PerLeave : public GeoObject {
public:
    //! constructor
    GeoObject_PerLeave() = delete;
    //! constructor
    GeoObject_PerLeave(TypeObject name_, Identifier identifier_, const vector<Identifier>& leaves_) : GeoObject(name_, identifier_, leaves_), hasPeriodicRoot{ false }, periodicRoot{ 0 } {}
    //! return whether the object is periodic or not
    bool isPeriodic() const { return hasPeriodicRoot; }
    //! return the identifier for sorting (when merging)
    IdentifierSort getId_forSort() const override;
    //! getter
    Identifier getPeriodicRoot() const { if (hasPeriodicRoot) return periodicRoot; else Merope_error_impossible(); }
    //! setter
    void setPeriodicRoot(const Identifier& periodicRoot_) { periodicRoot = periodicRoot_; hasPeriodicRoot = true; }
    //! setter
    void removePeriodicRoot() { hasPeriodicRoot = false; }
    //! shift the indices by a given shift
    void shiftIndices(const Identifier shift) override;
private:
    //! indicates, if applicable, if the object is periodic
    bool hasPeriodicRoot;
    //! links the object to is periodic father
    Identifier periodicRoot;
};

template<unsigned short DIM>
class GeoPoint final : public GeoObject_PerLeave {
public:
    Point<DIM> coordinates;
    //! constructor
    GeoPoint(Identifier id_, const Point<DIM>& coordinates_) : GeoObject_PerLeave(TypeObject::Point, id_, {}), coordinates(coordinates_) {}
    //! are two points the same?
    AreSame areSame(const GeoPoint<DIM>& point2, double epsilon) const {
        return static_cast<AreSame>(geomTools::distanceCarre<DIM>(this->coordinates, point2.coordinates) < epsilon * epsilon);
    }
    //! are two points the same in periodic coordinates
    AreSame areSamePer(const AmbiantSpace::Tore<DIM>& torus, const GeoPoint<DIM>& point2, double epsilon) const {
        return static_cast<AreSame>(torus.distanceCarre(this->coordinates, point2.coordinates) < epsilon * epsilon);
    }
    //! print the thing
    void print(std::ostream& f) const override { GeoObject_PerLeave::print(f); auxi_function::writeVectorToString(this->coordinates, f); f << endl; }
    //! a point is never singular
    bool isSingular() const override { return false; }
};


enum class TypeEdge {
    Segment, Circle
};

class Edge final : public GeoObject {
public:
    //! constructor
    Edge(Identifier identifier_, const vector<Identifier>& leaves_, TypeEdge typeEdge_ = TypeEdge::Segment)
        : GeoObject(TypeObject::Edge, identifier_, leaves_), typeEdge{ typeEdge_ } {}
    //! return whether the line is singular or not
    bool isSingular() const override {
        if (typeEdge == TypeEdge::Segment) return leaves[0] == leaves[1];
        Merope_error_impossible();
    }
    //! return whether the line is coherent or not
    bool isCoherent() const override {
        if (typeEdge == TypeEdge::Segment) return leaves.size() == 2;
        if (typeEdge == TypeEdge::Circle) return leaves.size() == 3;
        Merope_error_impossible();
    }
    TypeEdge typeEdge;
};


class CurveLoop final : public GeoObject {
public:
    //! constructor
    CurveLoop(Identifier identifier_, const vector<Identifier>& leaves_) : GeoObject(TypeObject::CurveLoop, identifier_, leaves_) {}
    //! test whether the object is singular or not
    bool isSingular() const override {
        bugTest();
        return (leaves.size() == 0) or (leaves.size() == 2 and leaves[0] == -leaves[1]);
    }
};

enum class TypeSurface {
    Plane, Curved
};

class Surface final : public GeoObject_PerLeave {
public:
    //! constructor
    Surface(Identifier identifier_, const vector<Identifier>& leaves_, TypeSurface typeSurface_ = TypeSurface::Plane) :
        GeoObject_PerLeave(TypeObject::Surface, identifier_, leaves_),
        typeSurface{ typeSurface_ } {}
    //! type of surface
    TypeSurface typeSurface;
};


class SurfaceLoop final : public GeoObject {
public:
    //! constructor
    SurfaceLoop(Identifier identifier_, const vector<Identifier>& leaves_) : GeoObject(TypeObject::SurfaceLoop, identifier_, leaves_) {}
};


class Solid final : public GeoObject {
public:
    //! constructor
    Solid(Identifier identifier_, const vector<Identifier>& leaves_) : GeoObject(TypeObject::Solid, identifier_, leaves_) {}
};

class PhysicalSurface final : public GeoObject {
public:
    //! constructor
    PhysicalSurface(Identifier identifier_, const vector<Identifier>& leaves_) : GeoObject(TypeObject::PhysicalSurface, identifier_, leaves_) {}
};

class PhysicalVolume final : public GeoObject {
public:
    //! constructor
    PhysicalVolume(Identifier identifier_, const vector<Identifier>& leaves_) : GeoObject(TypeObject::PhysicalSolid, identifier_, leaves_) {}
};

class PerPoint : public PerGeoObject {
public:
    //! constructor
    PerPoint(Identifier identifier_, const vector<Identifier>& leaves_) : PerGeoObject(TypeObject::PerPoint, identifier_, leaves_) {}
};


template<unsigned short DIM>
class PerSurface final : public PerGeoObject {
public:
    //! constructor
    //PerSurface(Identifier identifier_, const vector<Identifier>& leaves_): PerGeoObject(TypeObject::PeriodicSurface, identifier_, leaves_){}
    //! constructor
    PerSurface(Identifier identifier_, const vector<Identifier>& leaves_, const Point<DIM>& translation_) :
        PerGeoObject(TypeObject::PeriodicSurface, identifier_, leaves_), translation(translation_) {}
    //! translation between the first and second surface (=pt_surf_1 - pt_surf_2)
    Point<DIM> translation;
public:
    //! verify coherence
    bool isCoherent() const override { return leaves.size() == 2; }
    //! test whether the object is singular or not
    virtual bool isSingular()const override { return leaves.size() == 0; }
    //! swap the periodic surfaces
    void swapSurf();
};

//! remove all the singular entities
//! fixme
template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES>
bool findSingularAndRemove(DICT_ROOTS& dictRoots, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves);
//! remove all the singular entities
//! fixme
template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES, class DICT_PER_THINGS>
bool findSingularAndRemovePeriodic(DICT_ROOTS& dictRoots, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves, DICT_PER_THINGS& dictPerThings);
//! fixme
template<class DICT_THINGS>
vector<Identifier> findSingularThings(const DICT_THINGS& dictThings);
//! fixme
template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES>
void removeSingularThings(const vector<Identifier>& singularThings, DICT_ROOTS& dictRoots, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves);
//! fixme
template<class DICT_PER_THINGS, class DICT_THINGS>
void removePeriodicRoots(const vector<Identifier>& singularThings, DICT_PER_THINGS& dictPerThings, const DICT_THINGS& dictThings);
//! removes an object from the dictThings, and inside the dictLeaves
template<class DICT_THINGS, class DICT_LEAVES>
void removeObjectDownwards(Identifier idObject, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves);

//! test whether two geometric entities are the same or not
template<class OBJ>
AreSame areSame(const OBJ& obj1, const OBJ& obj2);

//! find all the singular shapes into a DICT_THINGS
template<class DICT_THINGS>
vector<typename DICT_THINGS::mapped_type> findSingular(DICT_THINGS dictThings);

//! tells whether the graph structure is coherent (leaves and roots are aware of themselves)
template<class DICT_THINGS, class DICT_LEAVES>
bool isGraphCoherent(const DICT_THINGS& dictThing, const DICT_LEAVES& dictLeaves);
//!
template<mesh::geoObjects::TypeRelation RELAT, class DICT_THINGS, class DICT_RELATED>
bool isGraphCoherent_auxi(const DICT_THINGS& dictThing, const DICT_RELATED& dictRelated);

//! merge list of two things
template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES>
void vec_merge(vector<SameThings<Identifier>> vecSameThings, DICT_ROOTS& dictRoots, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves);
//! merge two things
template<class DICT_ROOTS, class DICT_THINGS, class DICT_LEAVES>
void merge(const SameThings<Identifier> sameThings, DICT_ROOTS& dictRoots, DICT_THINGS& dictThings, DICT_LEAVES& dictLeaves);

}  // namespace geoObjects


namespace sameThings {
//! \see sameThings::getReplacementList here, implicitly use the function areSame
template<class DICT_THINGS>
vector<SameThings<Identifier>> getReplacementList(const DICT_THINGS& dictThing);
//! \see sameThings::getReplacementList here, implicitly use the function areSame
template<class DICT_THINGS, class COMPARISON_FUNCTION>
vector<SameThings<Identifier>> getReplacementList(const DICT_THINGS& dictThing, COMPARISON_FUNCTION comparisonFunction);
}  // namespace sameThings


}  // namespace mesh
}  // namespace merope

#include "../Mesh/GeoObjects.ixx"

