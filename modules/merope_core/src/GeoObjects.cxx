//! Copyright : see license.txt
//!
//! \briefEncodes geometry objects compatible with OpenCascad/gmsh


#include "Mesh/GeoObjects.hxx"
#include "../../GenericTools/AuxiFunctions.hxx"


namespace merope {
namespace mesh {
namespace geoObjects {

void GeoObject::removeLeafSingular(Identifier id) {
    auxi_function::erase_if(leaves, [&id](auto i) {return id == abs(i);});
}

void GeoObject::replaceLeaf(Identifier id_1, Identifier id_2, long sense) {
    for (auto& id : leaves) {
        if (abs(id) == abs(id_1)) {
            id = (id >= 0) ? id_2 * sense : -id_2 * sense;
        }
    }
}

sameThings::IdentifierSort GeoObject_PerLeave::getId_forSort() const {
    if (not isPeriodic()) {
        return IdentifierSort({ 0, identifier });
    } else {
        return IdentifierSort({ -periodicRoot, identifier });
    }
}

/*
void CurveLoop::simplify() {
    // fixme : diff from python, don't know why
    bool hasRemovedLeaf = false;
    for(size_t j = 0, k = 0; j < leaves.size(); j++){
        k = j+1 % leaves.size();
        if(leaves[j] == -leaves[k]){
            leaves.erase(leaves.begin() + max(j,k));
            leaves.erase(leaves.begin() + min(j,k));
            hasRemovedLeaf = true;
        }
    }
    if(hasRemovedLeaf){
        simplify();
    }
}
*/

void GeoObject::print(std::ostream& f) const {
    f << getName(name) << "(" << identifier << ") leaves = (";
    auxi_function::writeVectorToString(leaves, f);
    f << ")";
    f << " roots = (";
    auxi_function::writeVectorToString(roots, f);
    f << ")" << endl;
}

string getName(TypeObject typeObject) {
    static const vector<string> singleton{ "Unknown", "Point", "Edge", "CurveLoop", "Surface", "SurfaceLoop", "Solid", "PhysicalSurface", "PhysicalSolid", "PerPoint", "PeriodicSurface" };
    return singleton.at(static_cast<size_t>(typeObject));
}

void GeoObject::shiftIndices(const Identifier shift) {
    this->identifier += shift;
    //
    TYPE<TypeRelation::Root> newRoots{};
    for (auto& index : this->roots) {
        newRoots.insert(index + shift);
    }
    this->roots = newRoots;
    //
    for (auto& index : this->leaves) {
        if (index >= 0) index += shift;
        if (index < 0) index -= shift;
    }
}

void GeoObject_PerLeave::shiftIndices(const Identifier shift) {
    GeoObject::shiftIndices(shift);
    if (hasPeriodicRoot) {
        setPeriodicRoot(getPeriodicRoot() + shift);
    }
}

}  // namespace geoObjects
}  // namespace mesh
}  // namespace merope

