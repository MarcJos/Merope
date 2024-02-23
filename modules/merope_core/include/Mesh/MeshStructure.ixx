//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MESH_MESHSTRUCTURE_IXX_
#define MESH_MESHSTRUCTURE_IXX_


#include "../MeropeNamespace.hxx"


namespace merope {
namespace mesh {
namespace meshStructure {

// GeoStructure<DIM>

template<unsigned short DIM>
inline VoroMesh_NotPeriodic<DIM>::VoroMesh_NotPeriodic(VoroMesh_UnStructureData<DIM> rawData) :
    dictPoint{}, dictEdge{}, dictCurveLoop{}, dictSurface{}, dictSurfaceLoop{}, dictSolid{},
    torus(rawData.L) {
    translate(rawData.vecPoint, dictPoint);
    translate(rawData.vecEdge, dictEdge);
    translate(rawData.vecCurveLoop, dictCurveLoop);
    translate(rawData.vecSurface, dictSurface);
    translate(rawData.vecSurfaceLoop, dictSurfaceLoop);
    translate(rawData.vecSolid, dictSolid);
    this->buildTree();
}

template<unsigned short DIM>
inline void VoroMesh_NotPeriodic<DIM>::buildTree() {
    connectRoot(dictPoint, dictEdge);
    connectRoot(dictEdge, dictCurveLoop);
    connectRoot(dictCurveLoop, dictSurface);
    connectRoot(dictSurface, dictSurfaceLoop);
    connectRoot(dictSurfaceLoop, dictSolid);
}

template<unsigned short DIM>
inline bool VoroMesh_NotPeriodic<DIM>::isCoherent() const {
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
inline VoroMesh_Periodic<DIM>::VoroMesh_Periodic(VoroMesh_UnStructureData<DIM> rawData, double adim_epsilon_0, double adim_epsilon_1) :
    VoroMesh_NotPeriodic<DIM>(rawData), dictPerPoint{}, dictPerSurface{},
    epsilon_0{ adim_epsilon_0 * *(min_element(rawData.L.begin(), rawData.L.begin())) },
    epsilon_1{ adim_epsilon_1 * *(min_element(rawData.L.begin(), rawData.L.begin())) }{
    translate(rawData.vecPerSurface, dictPerSurface);
    unifyPoints();
    buildPeriodicity(epsilon_1);
    //! debug
    this->isCoherent();
    //!
    removeClosePoints();
    //! debug
    this->isCoherent();
}

template<unsigned short DIM>
inline bool VoroMesh_Periodic<DIM>::isStronglyCoherent() const {
    this->isCoherent();
    for (const auto& [id, surf] : this->dictSurface) {
        if ((surf.template get<TypeRelation::Root>().size() == 2 and not surf.isPeriodic())
            or (surf.template get<TypeRelation::Root>().size() == 1 and surf.isPeriodic())) {
            // OK
        } else {
            ofstream file("Errors.txt");
            this->print(file);
            cerr << "#############" << endl;
            cerr << "surface id" << id << endl;
            throw runtime_error("The structure has singular surfaces");
        }
    }
    return true;
}

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::restrictEnveloppe() {
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
    constexpr Identifier surfLoop_id = 1;
    this->dictSurfaceLoop = { {surfLoop_id, geoObjects::SurfaceLoop(surfLoop_id, allSurfaces)} };
    //!
    this->dictSolid = { {surfLoop_id, geoObjects::Solid(surfLoop_id, {surfLoop_id})} };
    //!
    this->buildTree();
    this->buildPeriodicity(this->epsilon_1);
}

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::unifyPoints() {
    auto listCouplePoints = tooClosePoints(this->epsilon_0);
    mergeAll(listCouplePoints);
    this->removeAllSingular();
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

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::buildPeriodicity(double epsilon_1_) {
    // reset periodic parts
    this->dictPerPoint = {};
    for (auto& [id, pt] : this->dictPoint) {
        pt.removePeriodicRoot();
    }
    //
    this->dictPerSurface = {};
    for (auto& [id, surf] : this->dictSurface) {
        surf.removePeriodicRoot();
    }
    //
    this->buildPerPoint(epsilon_1_);
    this->buildPerSurface();
    this->verifyPeriodicity();
}

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::verifyPeriodicity() {
    auto L = this->torus.L;
    auto goodOrdering = [&L](auto& pS) { // necessary to correct the mesh generation. Forces the first component of the translation to be negative
        size_t index = 0;
        for (; index < 3; index++) {
            if (abs(pS.translation[index]) > 0.5 * L[index]) {
                if (pS.translation[index] > 0) {
                    pS.swapSurf();
                }
                break;
            }
        }
        if (index == 3) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Unexpected");
        }
        };
    for (auto& [i, pS] : this->dictPerSurface) {
        auto ix = this->getPoints_from_Surface(pS.leaves[0]);
        auto iy = this->getPoints_from_Surface(pS.leaves[1]);
        if (ix.size() != iy.size()) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Unexpected : two periodic surfaces have different numbers of points");
        }
        goodOrdering(pS);
    }
}

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::removeClosePoints() {
    auto perTooClosePoints = this->periodicTooClosePoint(this->epsilon_1);
    auto& listCouplePoint = get<0>(perTooClosePoints);
    auto& listCouplePerPoint = get<1>(perTooClosePoints);
    // 1) if have to merge periodic points, deal with them first
    if (listCouplePerPoint.size() != 0) {
        this->mergePerPoints(listCouplePerPoint);
        this->removeClosePoints();
    }
    // 2) if not, deal with the others
    else if (listCouplePoint.size() != 0) {
        this->mergeAll(listCouplePoint);
    }
    this->removeAllSingular();
}

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::buildPerPoint(double epsilon_1_) {
    auto epsilon_0_copy = epsilon_1_;
    auto tore = this->torus;
    auto comparisonFunction = [epsilon_0_copy, &tore](const GeoPoint<DIM>& p1, const GeoPoint<DIM>& p2) {return p1.areSamePer(tore, p2, epsilon_0_copy);};
    //
    auto listPerPoint_id = sameThings::getReplacementList(this->dictPoint, comparisonFunction);
    for (auto& samePts : listPerPoint_id) {
        auto& pt1 = this->dictPoint.at(get<0>(samePts)), pt2 = this->dictPoint.at(get<1>(samePts));
        auto perPointIndentifier = pt2.identifier;
        if (dictPerPoint.find(perPointIndentifier) == dictPerPoint.end()) {
            dictPerPoint.insert(make_pair(perPointIndentifier, PerPoint(perPointIndentifier, { perPointIndentifier })));
        }
        dictPerPoint.at(perPointIndentifier).leaves.push_back(pt1.identifier);
        pt1.setPeriodicRoot(perPointIndentifier);
        pt2.setPeriodicRoot(perPointIndentifier);
    }
}

template<unsigned short DIM>
inline vector<Identifier> mesh::meshStructure::VoroMesh_Periodic<DIM>::getPoints_from_Surface(Identifier surf_id) const {
    std::vector<Identifier> result;
    auto cLoop_id = this->dictSurface.at(surf_id).leaves[0];
    for (auto edge_id : this->dictCurveLoop.at(cLoop_id).leaves) {
        if (edge_id > 0) result.push_back(this->dictEdge.at(abs(edge_id)).leaves[0]);
        else            result.push_back(this->dictEdge.at(abs(edge_id)).leaves[1]);
    }
    return result;
}

template<unsigned short DIM>
inline std::set<Identifier> VoroMesh_Periodic<DIM>::getSurfaces_from_Point(Identifier pt_id) const {
    std::set<Identifier> result{};
    for (Identifier edge_id : this->dictPoint.at(pt_id).template get<geoObjects::TypeRelation::Root>()) {
        for (Identifier cloop_id : this->dictEdge.at(edge_id).template get<geoObjects::TypeRelation::Root>()) {
            result.insert(*(this->dictCurveLoop.at(cloop_id).template get<geoObjects::TypeRelation::Root>().begin()));
        }
    }
    return result;
}

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::buildPerSurface() {
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
    for (auto& ps : dictPerSurface) {
        auto& perPlane = ps.second;
        for (auto idPs : perPlane.leaves) {
            this->dictSurface.at(idPs).setPeriodicRoot(perPlane.identifier);
        }
    }
}

template<unsigned short DIM>
inline Identifier VoroMesh_Periodic<DIM>::getMaxIndex() const {
    Identifier result = 0;
    auto getMax = [&result](const auto& dict) {
        if (dict.size() > 0) {
            auto elem = std::max_element(dict.begin(), dict.end(), [](const auto& obj1, const auto& obj2) {return obj1.first < obj2.first;});
            result = max(elem->first, result);
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
inline bool VoroMesh_Periodic<DIM>::comparePerSurface(Identifier surf_id1, Identifier surf_id2, const Point<DIM>& translation) const {
    // same surfaces are no perodic copies
    if (surf_id1 != surf_id2 and geomTools::normeCarre<DIM>(translation) > 0.5 * *(min_element(this->torus.L.begin(), this->torus.L.begin()))) {
        //
        vector<Identifier> pts_1 = getPoints_from_Surface(surf_id1);
        vector<Identifier> pts_2 = getPoints_from_Surface(surf_id2);
        size_t number_of_same_points = 0;
        for (auto it1 = pts_1.begin(); it1 != pts_1.end(); it1++) {
            for (auto it2 = pts_2.begin(); it2 != pts_2.end(); it2++) {
                auto& pt1 = this->dictPoint.at(*it1);
                auto& pt2 = this->dictPoint.at(*it2);
                if (geomTools::normeCarre<DIM>(translation - pt1.coordinates + pt2.coordinates) < auxi_function::puissance<2>(this->epsilon_0)) {
                    number_of_same_points++;
                    if (number_of_same_points > 2) {
                        return true;
                    }
                }
            }
        }
    }
    return false;
}

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::removeAllSingular() {
    bool hasRemoved = false;
    map<Identifier, GeoPoint<3>> emptyMap{};
    hasRemoved = hasRemoved or findSingularAndRemove(this->dictCurveLoop, this->dictEdge, this->dictPoint);
    hasRemoved = hasRemoved or findSingularAndRemove(this->dictSurface, this->dictCurveLoop, this->dictEdge);
    hasRemoved = hasRemoved or findSingularAndRemovePeriodic(this->dictSurfaceLoop, this->dictSurface, this->dictCurveLoop, this->dictPerSurface);
    hasRemoved = hasRemoved or findSingularAndRemove(this->dictSolid, this->dictSurfaceLoop, this->dictSurface);
    hasRemoved = hasRemoved or findSingularAndRemove(emptyMap, this->dictSolid, this->dictSurfaceLoop);
    if (hasRemoved) this->removeAllSingular();
}

template<unsigned short DIM>
inline vector<SameThings<Identifier>> mesh::meshStructure::VoroMesh_Periodic<DIM>::getClosePerPoints(
    const vector<SameThings<GeoPoint<DIM>>>& closePointDoublePer) const {
    auto dictPP = this->dictPerPoint;
    auto getPerPoints = [&dictPP](const auto& samePt) {
        auto id_perPoint1 = get<0>(samePt).getPeriodicRoot();
        auto id_perPoint2 = get<1>(samePt).getPeriodicRoot();
        return make_tuple(dictPP.at(id_perPoint1), dictPP.at(id_perPoint2), AreSame::Same);
        };
    vector<SameThings<PerPoint>> samePerPoints{};
    std::transform(closePointDoublePer.begin(), closePointDoublePer.end(), std::back_inserter(samePerPoints), getPerPoints);
    return sameThings::auxi::replaceGraph(samePerPoints); // removes duplicates and order well
}

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::mergePerPoints(const vector<SameThings<Identifier> >& samePoints) {
    for (const auto& sp : samePoints) {
        this->mergeSinglePerPoint(this->dictPerPoint.at(get<0>(sp)), this->dictPerPoint.at(get<1>(sp)));
    }
}

template<unsigned short DIM>
inline void mesh::meshStructure::VoroMesh_Periodic<DIM>::mergeSinglePerPoint(
    const PerPoint& perPt1, const PerPoint& perPt2) {
    auto refPoint = this->dictPoint.at(perPt1.leaves[0]);
    //! recover all the idPoints
    vector<Identifier>  allIdPoints = perPt1.leaves; allIdPoints.insert(allIdPoints.end(), perPt2.leaves.begin(), perPt2.leaves.end());
    sort(allIdPoints.begin(), allIdPoints.end());
    allIdPoints.erase(unique(allIdPoints.begin(), allIdPoints.end()), allIdPoints.end());
    //!
    for (const auto& idPt : allIdPoints) {
        auto& point = this->dictPoint.at(idPt);
        point.setPeriodicRoot(perPt2.identifier);
    }
    dictPerPoint.at(perPt2.identifier).leaves = allIdPoints;
    dictPerPoint.erase(perPt1.identifier);
    this->alignPerPoints(perPt2.identifier);
}

template<unsigned short DIM>
inline geoObjects::PhysicalSurface mesh::meshStructure::VoroMesh_Periodic<DIM>::getOuterSurface(Identifier id) const {
    vector<Identifier> allPerSurfaces{};
    for (const auto& [id_, perSurf] : this->dictPerSurface) {
        allPerSurfaces.push_back(perSurf.leaves[0]);
        allPerSurfaces.push_back(perSurf.leaves[1]);
    }
    sort(allPerSurfaces.begin(), allPerSurfaces.end());
    return geoObjects::PhysicalSurface(id, allPerSurfaces);
}

template<unsigned short DIM>
inline void mesh::meshStructure::VoroMesh_NotPeriodic<DIM>::print(std::ostream& f) const {
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
inline void mesh::meshStructure::VoroMesh_Periodic<DIM>::print(std::ostream& f) const {
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
inline void VoroMesh_Periodic<DIM>::alignPerPoints(Identifier perPoint_id) {
    const auto& list_pt_id = this->dictPerPoint.at(perPoint_id).leaves;
    if (list_pt_id.size() > 1) {
        const auto& refCoordinates = this->dictPoint.at(list_pt_id[0]).coordinates;
        Point<DIM>& L = this->torus.L;
        for (Identifier pt_id : list_pt_id) {
            auto& coordinates = this->dictPoint.at(pt_id).coordinates;
            for (size_t i = 0; i < DIM; i++) {
                coordinates[i] = round((coordinates[i] - refCoordinates[i]) / L[i]) * L[i] + refCoordinates[i];
            }
        }
    }
}

template<unsigned short DIM>
inline void VoroMesh_Periodic<DIM>::mergeAll(
    vector<SameThings<Identifier>> listSamePoints) {
    map<Identifier, GeoPoint<DIM>> emptyMap{};
    updatePeriodicMerge(listSamePoints, this->dictPoint, this->dictPerPoint);
    vec_merge(listSamePoints, this->dictEdge, this->dictPoint, emptyMap);
    //
    auto listSameEdges = sameThings::getReplacementList(this->dictEdge);
    vec_merge(listSameEdges, this->dictCurveLoop, this->dictEdge, this->dictPoint);
    //
    auto listSameCurveLoops = sameThings::getReplacementList(this->dictCurveLoop);
    vec_merge(listSameCurveLoops, this->dictSurface, this->dictCurveLoop, this->dictEdge);
    //
    auto listSameSurfaces = sameThings::getReplacementList(this->dictSurface);
    // ugly
    for (auto& sSurf : listSameSurfaces) {
        if (this->dictSurface.at(get<1>(sSurf)).leaves.at(0) < 0) { //  avoid negative surfaces
            swap(get<0>(sSurf), get<1>(sSurf));
        }
    }
    // ugly
    updatePeriodicMerge(listSameSurfaces, this->dictSurface, this->dictPerSurface);
    vec_merge(listSameSurfaces, this->dictSurfaceLoop, this->dictSurface, this->dictCurveLoop);
    //
    auto listSameSurfaceLoops = sameThings::getReplacementList(this->dictSurfaceLoop);
    vec_merge(listSameSurfaceLoops, this->dictSolid, this->dictSurfaceLoop, this->dictSurface);
}

template<unsigned short DIM>
inline vector<sameThings::SameThings<Identifier>> VoroMesh_Periodic<DIM>::tooClosePoints(double adimensionalDistance) const {
    double epsilon = adimensionalDistance * (*(min_element(this->torus.L.begin(), this->torus.L.end())));
    auto comparisonFunction = [epsilon](const GeoPoint<DIM>& p1, const GeoPoint<DIM>& p2) {return p1.areSame(p2, epsilon);};
    return sameThings::getReplacementList(this->dictPoint, comparisonFunction);
}

template<unsigned short DIM>
inline tuple<vector<sameThings::SameThings<Identifier>>, vector<sameThings::SameThings<Identifier>>>
VoroMesh_Periodic<DIM>::periodicTooClosePoint(double adimensionalDistance) const {
    const auto& dictPt = this->dictPoint;
    auto areBothPeriodic = [&dictPt](const auto& sameThing) {
        // if both are in the same periodic point, this is false
        const auto& pt1 = dictPt.at(get<0>(sameThing)), pt2 = dictPt.at(get<1>(sameThing));
        if (pt1.isPeriodic() and not pt2.isPeriodic()) {
            cerr << __PRETTY_FUNCTION__ << endl;
            throw runtime_error("Unexpected : two points are close but one is not periodic!");
        }
        return pt1.isPeriodic() and pt2.isPeriodic() and (pt1.getPeriodicRoot() != pt2.getPeriodicRoot());
        };
    auto sameConcretePoints = [&dictPt](const auto& sameThing) {
        return make_tuple(dictPt.at(get<0>(sameThing)), dictPt.at(get<1>(sameThing)), get<2>(sameThing));
        };
    //!
    auto closePoints = this->tooClosePoints(adimensionalDistance);
    vector<SameThings<GeoPoint<DIM>>> ptsNonPer{}, ptsDoublePer{};
    for (const auto& same_cp_id : closePoints) {
        if (areBothPeriodic(same_cp_id)) {
            ptsDoublePer.push_back(sameConcretePoints(same_cp_id));
        } else {
            ptsNonPer.push_back(sameConcretePoints(same_cp_id));
        }
    }
    return make_tuple(sameThings::auxi::replaceGraph(ptsNonPer), this->getClosePerPoints(ptsDoublePer));
}

// GeoUnStructureData<DIM>

template<unsigned short DIM>
inline Identifier VoroMesh_UnStructureData<DIM>::getMaxIndex() const {
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
inline void VoroMesh_UnStructureData<DIM>::shiftIndices(Identifier shift) {
    auto shiftLoc = [shift](auto& list) {
        for (auto& elem : list) {
            elem.shiftIndices(shift);
        }
        };
    applyOnAllVectors(shiftLoc);
}

template<unsigned short DIM>
inline void VoroMesh_UnStructureData<DIM>::reset() {
    applyOnAllVectors([](auto& vec) { vec = {};});
}

template<unsigned short DIM>
template<class FUNCTION>
inline void VoroMesh_UnStructureData<DIM>::applyOnAllVectors(FUNCTION function) {
    function(vecPoint);
    function(vecEdge);
    function(vecCurveLoop);
    function(vecSurface);
    function(vecSurfaceLoop);
    function(vecSolid);
    function(vecPerSurface);
}

template<unsigned short DIM>
template<class FUNCTION>
inline void VoroMesh_UnStructureData<DIM>::applyOnAllVectors(FUNCTION function) const {
    function(vecPoint);
    function(vecEdge);
    function(vecCurveLoop);
    function(vecSurface);
    function(vecSurfaceLoop);
    function(vecSolid);
    function(vecPerSurface);
}

template<unsigned short DIM>
inline void VoroMesh_UnStructureData<DIM>::fromGeoPerStructure(
    VoroMesh_Periodic<DIM> geoStructure) {
    this->L = geoStructure.torus.L;
    this->reset();
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
}

template<unsigned short DIM>
inline void VoroMesh_UnStructureData<DIM>::append_with_shift(VoroMesh_UnStructureData<DIM> rawData) {
    // fixme : unefficient
    auto index = this->getMaxIndex();
    // fixme : unefficient
    rawData.shiftIndices(index + 1);
    //
    this->append(rawData);
}

template<unsigned short DIM>
inline void VoroMesh_UnStructureData<DIM>::append(const VoroMesh_UnStructureData<DIM>& rawData) {
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
}

//! template functions

template<class DICT_THINGS, class DICT_PERIODIC_THINGS>
inline void updatePeriodicMerge(vector<SameThings<Identifier>> vecThings_id, const DICT_THINGS& dictThings, DICT_PERIODIC_THINGS& dictPerThings) {
    for (auto& sameThings : vecThings_id) {
        auto& thing1 = dictThings.at(get<0>(sameThings)), thing2 = dictThings.at(get<1>(sameThings));
        if (thing1.isPeriodic()) {
            auto& perThing = dictPerThings.at(thing1.getPeriodicRoot());
            perThing.replaceLeaf(thing1.identifier, thing2.identifier, 1);
        }
    }
}

template<class VEC, class DICT>
inline void translate(const VEC& vec, DICT& dict) {
    dict = {};
    set<Identifier> indexes{};
    for (const auto& elem : vec) {
        // verify coherence
        if (indexes.find(elem.identifier) != indexes.end()) {
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
        indexes.insert(elem.identifier);
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

template<unsigned short DIM>
inline bool verifyTranslate(const AmbiantSpace::Tore<DIM>& torus, const Point<DIM>& translation, double epsilon) {
    return geomTools::normeCarre(torus.projection(translation)) < epsilon * epsilon;
}

template<class DICT_LEAF, class DICT_ROOT>
void restrictTo_RootLeaves_withoutRootConnection(DICT_LEAF& dictThings, const DICT_ROOT& dictRoots) {
    // retrieve all the identifiers of things
    std::set<Identifier> all_things_id{};
    for (const auto& [id, root] : dictRoots) {
        for (auto thing_id : root.leaves) {
            all_things_id.insert(thing_id);
        }
    }
    auxi_function::erase_if(dictThings, [&all_things_id](const auto& thing) {
        return all_things_id.find(thing.second.identifier) == all_things_id.end();
        });
}


} // namespace meshStructure
} // namespace mesh
} // namespace merope

#endif /* MESH_MESHSTRUCTURE_IXX_ */
