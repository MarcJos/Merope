//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../GenericTools/CPP_Functions.hxx"

#include "../../../AlgoPacking/include/AlgoRSA.hxx"
#include "../../../AlgoPacking/include/MultiDArrayObject.hxx"


namespace merope {
namespace mesh {
namespace meshStructure {

// GeoStructure<DIM>

template<unsigned short DIM>
VoroMesh_NotPeriodic<DIM>::VoroMesh_NotPeriodic(const VoroMesh_UnStructureData<DIM>& rawData) :
    dictPoint{}, dictEdge{}, dictCurveLoop{}, dictSurface{}, dictSurfaceLoop{}, dictSolid{},
    torus(rawData.L) {
    this->add_raw_mesh_data(rawData);
}

template<unsigned short DIM>
void VoroMesh_NotPeriodic<DIM>::add_raw_mesh_data(const VoroMesh_UnStructureData<DIM>& rawData) {
    translate(rawData.vecPoint, dictPoint);
    translate(rawData.vecEdge, dictEdge);
    translate(rawData.vecCurveLoop, dictCurveLoop);
    translate(rawData.vecSurface, dictSurface);
    translate(rawData.vecSurfaceLoop, dictSurfaceLoop);
    translate(rawData.vecSolid, dictSolid);
    this->buildTree();
}

template<unsigned short DIM>
void VoroMesh_NotPeriodic<DIM>::buildTree() {
    connectRoot(dictPoint, dictEdge);
    connectRoot(dictEdge, dictCurveLoop);
    connectRoot(dictCurveLoop, dictSurface);
    connectRoot(dictSurface, dictSurfaceLoop);
    connectRoot(dictSurfaceLoop, dictSolid);
}

template<unsigned short DIM>
bool VoroMesh_NotPeriodic<DIM>::isCoherent() const {
    auto testCoherence_loc = [this](const auto& dictRoot, const auto& dictLeaves, string nameError) {
        if (not geoObjects::isGraphCoherent(dictRoot, dictLeaves)) {
            this->print(cerr);
            string errorDisplay = string("Coherence Problem : ") + nameError;
            throw runtime_error(errorDisplay);
        }
        };
    testCoherence_loc(dictSolid, dictSurfaceLoop, "DictSolid");
    testCoherence_loc(dictSurfaceLoop, dictSurface, "DictSurfaceLoop");
    testCoherence_loc(dictSurface, dictCurveLoop, "DictSurface");
    testCoherence_loc(dictCurveLoop, dictEdge, "DictCurveLoop");
    testCoherence_loc(dictEdge, dictPoint, "DictEdge");
    return true;
}


// GeoPerStructure<DIM>

template<unsigned short DIM>
VoroMesh_Periodic<DIM>::VoroMesh_Periodic(const VoroMesh_UnStructureData<DIM>& rawData, bool check_coherence) :
    VoroMesh_NotPeriodic<DIM>(rawData), dictPerPoint{}, dictPerSurface{}{
    translate(rawData.vecPerPoint, dictPerPoint);
    translate(rawData.vecPerSurface, dictPerSurface);
    connectPerRoot(this->dictPoint, dictPerPoint);
    connectPerRoot(this->dictSurface, dictPerSurface);
    //
    restrictTo_RootLeaves_withoutRootConnection(this->dictCurveLoop, this->dictSurface);
    restrictTo_RootLeaves_withoutRootConnection(this->dictEdge, this->dictCurveLoop);
    restrictTo_RootLeaves_withoutRootConnection(this->dictPoint, this->dictEdge);
    //
    bool hasChangedStructure = true;
    while (hasChangedStructure) {
        hasChangedStructure = false;
        hasChangedStructure = this->removeAllSingular() or hasChangedStructure;
        hasChangedStructure = this->mergeAll() or hasChangedStructure;
    }
    this->buildPerSurface();
    this->orderPerSurface();
    // verifications
    if (check_coherence) {
        this->verifyPeriodicity();
        this->isStronglyCoherent();
    }
}

template<unsigned short DIM>
bool VoroMesh_Periodic<DIM>::isStronglyCoherent() const {
    // test local coherence
    auto local_coherence = [](const auto& dictThing) {
        for (const auto& [id, thing] : dictThing) {
            if (not thing.isCoherent()) {
                std::cerr << __PRETTY_FUNCTION__ << endl;
                thing.print(std::cerr); std::cerr << endl;
                throw runtime_error("Incoherent state");
            }
        }
        };
    local_coherence(this->dictPoint);
    local_coherence(this->dictEdge);
    local_coherence(this->dictCurveLoop);
    local_coherence(this->dictSurface);
    local_coherence(this->dictSurfaceLoop);
    local_coherence(this->dictSolid);
    local_coherence(this->dictPerPoint);
    local_coherence(this->dictPerSurface);

    this->isCoherent();

    // test if surfaces are shared
    for (const auto& [id, surf] : this->dictSurface) {
        if ((surf.template get<TypeRelation::Root>().size() == 2 and not surf.isPeriodic())
            or (surf.template get<TypeRelation::Root>().size() == 1 and surf.isPeriodic())
            or (surf.template get<TypeRelation::Root>().size() == 1 and this->dictSurfaceLoop.at(*(surf.template get<TypeRelation::Root>().begin())).template get<TypeRelation::Root>().size() == 2)) {
            // OK
        } else {
            ofstream file("Errors.txt");
            this->print(file);
            cerr << "#############" << endl;
            cerr << "surface id : " << id << endl;
            cerr << "number of parents : " << surf.template get<TypeRelation::Root>().size() << endl;
            cerr << std::boolalpha << "is periodic : " << surf.isPeriodic() << endl;
            cerr << "see Errors.txt for infos on the mesh" << endl;
            throw runtime_error("The structure has singular surfaces");
        }
    }
    return true;
}

template<unsigned short DIM>
void VoroMesh_Periodic<DIM>::restrictEnveloppe() {
    restrictTo_RootLeaves_withoutRootConnection(this->dictSurface, this->dictPerSurface);
    restrictTo_RootLeaves_withoutRootConnection(this->dictCurveLoop, this->dictSurface);
    restrictTo_RootLeaves_withoutRootConnection(this->dictEdge, this->dictCurveLoop);
    restrictTo_RootLeaves_withoutRootConnection(this->dictPoint, this->dictEdge);
    //!
    // Build the whole surfaceLoop
    vector<Identifier> allSurfaces{};
    for (const auto& [id, surf] : this->dictSurface) {
        allSurfaces.push_back(id);
    }
    //!
    constexpr Identifier surfLoop_id = 1;
    this->dictSurfaceLoop = { {surfLoop_id, geoObjects::SurfaceLoop(surfLoop_id, allSurfaces)} };
    this->dictSolid = { {surfLoop_id, geoObjects::Solid(surfLoop_id, {surfLoop_id})} };
    //!
    this->buildTree();
}



template<unsigned short DIM>
void VoroMesh_Periodic<DIM>::orderPerSurface() {
    auto L = this->torus.L;
    auto goodOrdering = [&L](auto& pS) {  // necessary to correct the mesh generation. Forces the first component of the translation to be negative
        for (size_t index = 0; index < 3; index++) {
            if (abs(pS.translation[index]) > 0.5 * L[index]) {
                if (pS.translation[index] > 0) {
                    pS.swapSurf();
                }
                break;
            }
        }
        };
    for (auto& [i, pS] : this->dictPerSurface) {
        goodOrdering(pS);
    }
}

template<unsigned short DIM>
void VoroMesh_Periodic<DIM>::verifyPeriodicity() const {
    for (auto& [i, pS] : this->dictPerSurface) {
        auto ix = this->getPoints_from_Surface(pS.leaves[0]);
        auto iy = this->getPoints_from_Surface(pS.leaves[1]);
        Merope_assert(ix.size() == iy.size(),
            "Unexpected : two periodic surfaces have different numbers of points");
        // check if all points coincide one to one with periodicity
        bool correct_surface = true;
        auto display_error = [&](int number) {
            cerr << __PRETTY_FUNCTION__ << endl;
            cerr << "Surface " << i << endl;
            cerr << "warning : " << number << endl;
            cerr << "All points of periodic surfaces are not identified as periodic points." << endl;
            cerr << "Attempting to correct it" << endl;
            correct_surface = false;
            };
        for (auto it1 = ix.begin(); it1 != ix.end(); it1++) {
            auto& pt1 = this->dictPoint.at(*it1);
            if (not pt1.isPeriodic()) {
                display_error(1);
                break;
            }
            //
            bool success = false;
            for (auto it2 = iy.begin(); it2 != iy.end(); it2++) {
                auto& pt2 = this->dictPoint.at(*it2);
                if (not pt2.isPeriodic()) {
                    display_error(2);
                    break;
                } else if (pt1.getPeriodicRoot() == pt2.getPeriodicRoot()) {
                    success = true;
                    break;
                }
            }
            //
            if (not success) {
                display_error(3);
            }
        }
        // if not all points coincide one to one with periodicity, try to correct
        if (not correct_surface) {
            bool success = false;
            long index_s_1, index_s_2;
            for (index_s_1 = 0; index_s_1 != ix.size(); index_s_1++) {
                auto& pt1 = this->dictPoint.at(ix[index_s_1]);
                if (pt1.isPeriodic()) {
                    break;
                }
            }
            for (index_s_2 = 0; index_s_2 != iy.size(); index_s_2++) {
                auto& pt1 = this->dictPoint.at(ix[index_s_1]);
                auto& pt2 = this->dictPoint.at(iy[index_s_2]);
                if (pt2.isPeriodic() and pt1.getPeriodicRoot() == pt2.getPeriodicRoot()) {
                    success = true;
                    break;
                }
            }
            if (not success) {
                cerr << __PRETTY_FUNCTION__ << endl;
                throw runtime_error("Unexpected that periodic surfaces share no common periodic point");
            }
            for (int j = 0; j < ix.size(); j++) {
                index_s_1++;
                index_s_2--;
                index_s_1 = auxi_function::fast_modulo(index_s_1, ix.size());
                index_s_2 = auxi_function::fast_modulo(index_s_2, ix.size());
                auto& pt1 = this->dictPoint.at(ix[index_s_1]);
                auto& pt2 = this->dictPoint.at(iy[index_s_2]);
                Merope_assert((pt1.isPeriodic() and pt2.isPeriodic()) and (pt1.getPeriodicRoot() != pt2.getPeriodicRoot()), "Two points are not recognized as periodic");
            }
        }
    }
}

template<unsigned short DIM>
vector<Identifier> mesh::meshStructure::VoroMesh_Periodic<DIM>::getPoints_from_Surface(Identifier surf_id) const {
    std::vector<Identifier> result;
    auto cLoop_id = this->dictSurface.at(surf_id).leaves[0];
    for (auto edge_id : this->dictCurveLoop.at(cLoop_id).leaves) {
        if (edge_id > 0) result.push_back(this->dictEdge.at(abs(edge_id)).leaves[0]);
        else             result.push_back(this->dictEdge.at(abs(edge_id)).leaves[1]);
    }
    return result;
}

template<unsigned short DIM>
std::set<Identifier> VoroMesh_Periodic<DIM>::getSurfaces_from_Point(Identifier pt_id) const {
    std::set<Identifier> result{};
    for (Identifier edge_id : this->dictPoint.at(pt_id).template get<geoObjects::TypeRelation::Root>()) {
        for (Identifier cloop_id : this->dictEdge.at(edge_id).template get<geoObjects::TypeRelation::Root>()) {
            result.insert(*(this->dictCurveLoop.at(cloop_id).template get<geoObjects::TypeRelation::Root>().begin()));
        }
    }
    return result;
}

template<unsigned short DIM>
void VoroMesh_Periodic<DIM>::buildPerSurface() {
    //! get all candidate surfaces
    for (const auto& pp : dictPerPoint) {
        const auto& perPoint = pp.second;
        for (size_t i = 0; i < perPoint.leaves.size(); i++) {
            for (size_t j = i + 1; j < perPoint.leaves.size(); j++) {
                Identifier pt_id_1 = perPoint.leaves[i], pt_id_2 = perPoint.leaves[j];
                Point<DIM> translation = this->dictPoint.at(pt_id_1).coordinates - this->dictPoint.at(pt_id_2).coordinates;
                if (geomTools::normeCarre<DIM>(translation) > 0.5 * auxi_function::puissance<2>(*min_element(this->torus.L.begin(), this->torus.L.end()))) {
                    auto set_surf_id1 = this->getSurfaces_from_Point(pt_id_1);
                    auto set_surf_id2 = this->getSurfaces_from_Point(pt_id_2);
                    for (auto surf_id1 : set_surf_id1) {
                        for (auto surf_id2 : set_surf_id2) {
                            if (this->dictPerSurface.find(surf_id1) == this->dictPerSurface.end() and this->dictPerSurface.find(surf_id2) == this->dictPerSurface.end()) {
                                if (this->comparePerSurface(surf_id1, surf_id2, translation)) {
                                    Identifier id_perSurface = min(surf_id1, surf_id2);
                                    this->dictPerSurface.insert(make_pair(id_perSurface, PerSurface<DIM>(id_perSurface, { surf_id1, surf_id2 }, translation)));
                                }
                            }
                        }
                    }
                }
            }
        }
    }
    //! update dictSurface
    connectPerRoot(this->dictSurface, dictPerSurface);
}

template<unsigned short DIM>
Identifier VoroMesh_Periodic<DIM>::getMaxIndex() const {
    Identifier result = 0;
    auto getMax = [&result](const auto& dict) {
        if (dict.size() > 0) {
            auto elem = std::max_element(dict.begin(), dict.end(), [](const auto& obj1, const auto& obj2) {return obj1.first < obj2.first;});
            if (elem != dict.end()) result = max(elem->first, result);
        }
        };
    getMax(this->dictPoint);
    getMax(this->dictEdge);
    getMax(this->dictCurveLoop);
    getMax(this->dictSurface);
    getMax(this->dictSurfaceLoop);
    getMax(this->dictSolid);
    return result;
}

template<unsigned short DIM>
bool VoroMesh_Periodic<DIM>::comparePerSurface(Identifier surf_id1, Identifier surf_id2, const Point<DIM>& translation) const {
    // same surfaces are no perodic copies
    if (surf_id1 != surf_id2 and geomTools::normeCarre<DIM>(translation) > 0.5 * *(min_element(this->torus.L.begin(), this->torus.L.begin()))) {
        //
        vector<Identifier> pts_1 = getPoints_from_Surface(surf_id1);
        vector<Identifier> pts_2 = getPoints_from_Surface(surf_id2);
        if (pts_1.size() != pts_2.size()) {
            return false;
        }
        // magical ( warning )
        constexpr double tol = 1e-10;
        for (long i_begin = 0; i_begin < pts_1.size(); i_begin++) {
            bool success = true;
            for (long i = 0; i < pts_1.size(); i++) {
                long j = auxi_function::fast_modulo(i_begin - i, pts_1.size());
                if (geomTools::normeCarre<DIM>(this->dictPoint.at(pts_1[i]).coordinates - this->dictPoint.at(pts_2[j]).coordinates - translation) > tol * tol * geomTools::normeCarre<DIM>(this->torus.L)) {
                    success = false;
                    break;
                }
            }
            if (success) {
                return true;
            }
        }
    }
    return false;
}

template<unsigned short DIM>
bool VoroMesh_Periodic<DIM>::removeAllSingular() {
    bool hasRemoved = false;
    map<Identifier, GeoPoint<3>> emptyMap{};
    hasRemoved = hasRemoved or findSingularAndRemove(this->dictCurveLoop, this->dictEdge, this->dictPoint);
    hasRemoved = hasRemoved or findSingularAndRemove(this->dictSurface, this->dictCurveLoop, this->dictEdge);
    hasRemoved = hasRemoved or findSingularAndRemovePeriodic(this->dictSurfaceLoop, this->dictSurface, this->dictCurveLoop, this->dictPerSurface);
    hasRemoved = hasRemoved or findSingularAndRemove(this->dictSolid, this->dictSurfaceLoop, this->dictSurface);
    hasRemoved = hasRemoved or findSingularAndRemove(emptyMap, this->dictSolid, this->dictSurfaceLoop);
    if (hasRemoved) {
        this->removeAllSingular();
    } else {
        // verifications
        // here, only close points related to opposed faces should be merged, not the ones with small edges
        auto singularEdges = geoObjects::findSingular(this->dictEdge);
        if (singularEdges.size() != 0) {
            for (const auto& edge : singularEdges) {
                cerr << edge.leaves[0] << endl;
                cerr << edge.identifier << endl;
            }
            throw runtime_error("Singular line");
        }
    }
    return hasRemoved;
}

template<unsigned short DIM>
geoObjects::PhysicalSurface mesh::meshStructure::VoroMesh_Periodic<DIM>::getPeriodicOuterSurface(Identifier id) const {
    vector<Identifier> allPerSurfaces{};
    for (const auto& [id_, perSurf] : this->dictPerSurface) {
        allPerSurfaces.push_back(perSurf.leaves[0]);
        allPerSurfaces.push_back(perSurf.leaves[1]);
    }
    sort(allPerSurfaces.begin(), allPerSurfaces.end());
    return geoObjects::PhysicalSurface(id, allPerSurfaces);
}

template<unsigned short DIM>
void mesh::meshStructure::VoroMesh_NotPeriodic<DIM>::print(std::ostream& f) const {
    f << "--------------" << endl;
    f << "GeoStructure<" << DIM << ">" << endl;
    auto printAll = [&f](const auto& dict) {
        for (const auto& obj : dict) {
            obj.second.print(f);
        }
        f << endl << endl;
        };
    printAll(this->dictPoint);
    printAll(this->dictEdge);
    printAll(this->dictCurveLoop);
    printAll(this->dictSurface);
    printAll(this->dictSurfaceLoop);
    printAll(this->dictSolid);
    f << "--------------" << endl;
}

template<unsigned short DIM>
void mesh::meshStructure::VoroMesh_Periodic<DIM>::print(std::ostream& f) const {
    f << "--------------" << endl;
    f << "Nb Point" << this->dictPoint.size() << endl;
    f << "Nb Edges" << this->dictEdge.size() << endl;
    f << "Nb CurveLoop" << this->dictCurveLoop.size() << endl;
    f << "Nb Surface" << this->dictSurface.size() << endl;
    f << "Nb SurfaceLoop" << this->dictSurfaceLoop.size() << endl;
    f << "Nb Solids" << this->dictSolid.size() << endl;
    f << "Nb Perpoints" << this->dictPerPoint.size() << endl;
    f << "Nb PerSurf" << this->dictPerSurface.size() << endl;
    f << "--------------" << endl;
    this->VoroMesh_NotPeriodic<DIM>::print(f);
    f << "--------------" << endl;
    auto printAll = [&f](const auto& dict) {
        for (const auto& obj : dict) {
            obj.second.print(f);
        }
        f << endl << endl;
        };
    printAll(this->dictPerPoint);
    printAll(this->dictPerSurface);
}

template<unsigned short DIM>
bool VoroMesh_Periodic<DIM>::mergeAll() {
    bool hasMerged = false;
    auto listSameEdges = sameThings::getReplacementList(this->dictEdge);
    vec_merge(listSameEdges, this->dictCurveLoop, this->dictEdge, this->dictPoint);
    hasMerged = hasMerged or (listSameEdges.size() > 0);
    //
    auto listSameCurveLoops = sameThings::getReplacementList(this->dictCurveLoop);
    vec_merge(listSameCurveLoops, this->dictSurface, this->dictCurveLoop, this->dictEdge);
    hasMerged = hasMerged or (listSameCurveLoops.size() > 0);
    //
    auto listSameSurfaces = sameThings::getReplacementList(this->dictSurface);
    // ugly
    for (auto& sSurf : listSameSurfaces) {
        if (this->dictSurface.at(get<1>(sSurf)).leaves.at(0) < 0) {  //  avoid negative surfaces
            swap(get<0>(sSurf), get<1>(sSurf));
        }
    }
    // ugly
    updatePeriodicMerge(listSameSurfaces, this->dictSurface, this->dictPerSurface);
    vec_merge(listSameSurfaces, this->dictSurfaceLoop, this->dictSurface, this->dictCurveLoop);
    hasMerged = hasMerged or (listSameSurfaces.size() > 0);
    //
    auto listSameSurfaceLoops = sameThings::getReplacementList(this->dictSurfaceLoop);
    vec_merge(listSameSurfaceLoops, this->dictSolid, this->dictSurfaceLoop, this->dictSurface);
    hasMerged = hasMerged or (listSameSurfaceLoops.size() > 0);
    return hasMerged;
}

// GeoUnStructureData<DIM>

template<unsigned short DIM>
Identifier VoroMesh_UnStructureData<DIM>::getMaxIndex() const {
    Identifier maxIndex = 0;
    auto getMax = [&maxIndex](const auto& list) {
        if (list.size() > 0) {
            auto elem = min_element(list.begin(), list.end(), [](const auto& e1, const auto& e2) {return e1.identifier > e2.identifier;});
            maxIndex = max(elem->identifier, maxIndex);
        }
        };
    applyOnAllVectors(getMax);
    return maxIndex;
}

template<unsigned short DIM>
void VoroMesh_UnStructureData<DIM>::shiftIndices(Identifier shift) {
    auto shiftLoc = [shift](auto& list) {
        for (auto& elem : list) {
            elem.shiftIndices(shift);
        }
        };
    applyOnAllVectors(shiftLoc);
}

template<unsigned short DIM>
void VoroMesh_UnStructureData<DIM>::reset() {
    applyOnAllVectors([](auto& vec) { vec = {};});
}

template<unsigned short DIM>
template<class FUNCTION>
void VoroMesh_UnStructureData<DIM>::applyOnAllVectors(FUNCTION function) {
    function(vecPoint);
    function(vecEdge);
    function(vecCurveLoop);
    function(vecSurface);
    function(vecSurfaceLoop);
    function(vecSolid);
    function(vecPerSurface);
    function(vecPerPoint);
}

template<unsigned short DIM>
template<class FUNCTION>
void VoroMesh_UnStructureData<DIM>::applyOnAllVectors(FUNCTION function) const {
    function(vecPoint);
    function(vecEdge);
    function(vecCurveLoop);
    function(vecSurface);
    function(vecSurfaceLoop);
    function(vecSolid);
    function(vecPerSurface);
    function(vecPerPoint);
}

template<unsigned short DIM>
void VoroMesh_UnStructureData<DIM>::fromGeoPerStructure(
    const VoroMesh_Periodic<DIM>& geoStructure) {
    this->reset();
    this->L = geoStructure.torus.L;
    auto copyInto = [](auto list1, auto& list2) {
        for (const auto& item : list1) {
            list2.push_back(item.second);
        }
        };
    copyInto(geoStructure.dictPoint, this->vecPoint);
    copyInto(geoStructure.dictEdge, this->vecEdge);
    copyInto(geoStructure.dictCurveLoop, this->vecCurveLoop);
    copyInto(geoStructure.dictSurface, this->vecSurface);
    copyInto(geoStructure.dictSurfaceLoop, this->vecSurfaceLoop);
    copyInto(geoStructure.dictSolid, this->vecSolid);
    copyInto(geoStructure.dictPerSurface, this->vecPerSurface);
    copyInto(geoStructure.dictPerPoint, this->vecPerPoint);
}

template<unsigned short DIM>
void VoroMesh_UnStructureData<DIM>::append_with_shift(VoroMesh_UnStructureData<DIM> rawData) {
    // fixme : unefficient
    auto index = this->getMaxIndex();
    // fixme : unefficient
    rawData.shiftIndices(index + 1);
    //
    this->append(rawData);
}

template<unsigned short DIM>
void VoroMesh_UnStructureData<DIM>::append(const VoroMesh_UnStructureData<DIM>& rawData) {
    auto append_aux = [](auto& vec1, const auto& vec2) {
        vec1.insert(vec1.end(), vec2.begin(), vec2.end());
        };
    append_aux(this->vecPoint, rawData.vecPoint);
    append_aux(this->vecEdge, rawData.vecEdge);
    append_aux(this->vecCurveLoop, rawData.vecCurveLoop);
    append_aux(this->vecSurface, rawData.vecSurface);
    append_aux(this->vecSurfaceLoop, rawData.vecSurfaceLoop);
    append_aux(this->vecSolid, rawData.vecSolid);
    append_aux(this->vecPerSurface, rawData.vecPerSurface);
    append_aux(this->vecPerPoint, rawData.vecPerPoint);
}

//! template functions

template<class DICT_THINGS, class DICT_PERIODIC_THINGS>
void updatePeriodicMerge(vector<SameThings<Identifier>> vecThings_id, const DICT_THINGS& dictThings, DICT_PERIODIC_THINGS& dictPerThings) {
    for (auto& sameThings : vecThings_id) {
        auto& thing1 = dictThings.at(get<0>(sameThings)), thing2 = dictThings.at(get<1>(sameThings));
        if (thing1.isPeriodic()) {
            auto& perThing = dictPerThings.at(thing1.getPeriodicRoot());
            perThing.replaceLeaf(thing1.identifier, thing2.identifier, 1);
        }
    }
}

template<class VEC, class DICT>
void translate(const VEC& vec, DICT& dict) {
    for (const auto& elem : vec) {
        // verify coherence
        if (dict.find(elem.identifier) != dict.end()) {
            cerr << __PRETTY_FUNCTION__ << endl;
            cerr << geoObjects::getName(elem.name) << endl;
            cerr << elem.identifier << endl;
            throw runtime_error("Cannot insert two elements with the same identifier");
        }
        if (elem.identifier <= 0) {
            cerr << __PRETTY_FUNCTION__ << endl;
            cerr << geoObjects::getName(elem.name) << endl;
            cerr << elem.identifier << endl;
            throw runtime_error("Cannot insert element with nonpositive identifier");
        }
        //
        dict.insert(make_pair(elem.identifier, elem));
    }
}

template<class DICT_THINGS, class DICT_ROOT>
void connectRoot(DICT_THINGS& dictThings, const DICT_ROOT& dictRoots) {
    for (auto& [id, leaf] : dictThings) {
        leaf.resetRoots();
    }
    for (const auto& [id, root] : dictRoots) {
        for (const auto& leaf : root.leaves) {
            dictThings.at(abs(leaf)).addRoot(id);
        }
    }
}

template<class DICT_THINGS, class DICT_ROOT>
void connectPerRoot(DICT_THINGS& dictThings, const DICT_ROOT& dictPerRoots) {
    for (auto& [id, leaf] : dictThings) {
        leaf.removePeriodicRoot();
    }
    for (auto& ps : dictPerRoots) {
        auto& perRoot = ps.second;
        for (auto idPs : perRoot.leaves) {
            dictThings.at(idPs).setPeriodicRoot(perRoot.identifier);
        }
    }
}

template<unsigned short DIM>
bool verifyTranslate(const AmbiantSpace::Tore<DIM>& torus, const Point<DIM>& translation, double epsilon) {
    return geomTools::normeCarre(torus.projection(translation)) < epsilon * epsilon;
}

template<class DICT_LEAF, class DICT_ROOT>
void restrictTo_RootLeaves_withoutRootConnection(DICT_LEAF& dictThings, const DICT_ROOT& dictRoots) {
    // retrieve all the identifiers of things
    std::set<Identifier> all_things_id{};
    for (const auto& [id, root] : dictRoots) {
        for (auto thing_id : root.leaves) {
            all_things_id.insert(abs(thing_id));
        }
    }
    auxi_function::erase_if(dictThings, [&all_things_id](const auto& thing) {
        return all_things_id.find(thing.second.identifier) == all_things_id.end();
        });
}


}  // namespace meshStructure
}  // namespace mesh
}  // namespace merope


