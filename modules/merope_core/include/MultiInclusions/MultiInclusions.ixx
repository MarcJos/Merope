//! Copyright : see license.txt
//!
//! \brief 
//!

#ifndef MULTIINCLUSIONS_IXX_
#define MULTIINCLUSIONS_IXX_

#include "../MeropeNamespace.hxx"


namespace merope {

template<unsigned short DIM>
template<class C>
inline vector<C>& MultiInclusions<DIM>::getInclusions() {
    if constexpr (is_same<C, smallShape::SphereInc<DIM>>::value) {
        return sphereInc;
    }
    if constexpr (is_same<C, smallShape::ConvexPolyhedronInc<DIM>>::value) {
        return polyhedrons;
    }
    if constexpr (is_same<C, smallShape::EllipseInc<DIM>>::value) {
        return ellipseInc;
    }
    if constexpr (is_same<C, smallShape::SpheroPolyhedronInc<DIM>>::value) {
        return spheroPolyhedrons;
    }
    cerr << __PRETTY_FUNCTION__ << endl;
    throw runtime_error("Unexpected");
}

template<unsigned short DIM>
template<class C>
inline void MultiInclusions<DIM>::setInnerShapes(const C& inclusions) {
    // reinitialize the shapes
    applyOnAllInclusions([](auto& vectorInclusions) {
        vectorInclusions = {};
        });
    // fill in the correct shape
    using INCLUSION_TYPE = typename C::value_type;
    vector<INCLUSION_TYPE>& vectorShapes = this->template getInclusions<INCLUSION_TYPE>();
    copy(inclusions.begin(), inclusions.end(), std::back_inserter(vectorShapes));
    sort(vectorShapes.begin(), vectorShapes.end(),
        [](const auto& poly1, const auto& poly2) {
            return poly1.identifier < poly2.identifier;
        });
}

template<unsigned short DIM>
inline MultiInclusions<DIM>::MultiInclusions() :
    InsideTorus<DIM>(), polyhedrons{ }, sphereInc{ }, ellipseInc{}, spheroPolyhedrons{}, matrixPresence{ false }, matrixPhase{ 0 } {
}

template<unsigned short DIM>
template<class C>
inline void MultiInclusions<DIM>::setInclusions_T(const C& vectorOfInclusions) {
    if constexpr (not(std::is_same<typename C::value_type, Sphere<DIM>>::value or std::is_same<typename C::value_type, Ellipse<DIM>>::value)) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected");
    }
    using INCLUSION_TYPE = typename std::conditional<std::is_same<typename C::value_type, Sphere<DIM>>::value, smallShape::SphereInc<DIM>, smallShape::EllipseInc<DIM>>::type;
    //!
    vector<INCLUSION_TYPE> theSphereInc{};
    transform(vectorOfInclusions.begin(), vectorOfInclusions.end(), std::back_inserter(theSphereInc), [](const auto& sph) {
        return INCLUSION_TYPE(sph);
        });
    for (size_t i = 0; i < theSphereInc.size(); i++) {
        theSphereInc[i].identifier = i;
    }
    setInnerShapes(theSphereInc);
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::setInclusions(LaguerreTess<DIM> polyX) {
    this->setLength(polyX.getL());
    polyX.computeTessels();
    setInnerShapes(polyX.getMicroInclusions());
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::setInclusions(const PolyInclusions<DIM>& polyInc) {
    this->setLength(polyInc.getL());
    this->setMatrixPhase(0);
    setInnerShapes(polyInc.getMicroInclusions());
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::setInclusions(
    const SphereInclusions<DIM>& sphereI) {
    this->setLength(sphereI.getL());
    this->setMatrixPhase(0);
    this->setInclusions_T(sphereI.getSpheres());
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::setInclusions(const EllipseInclusions<DIM>& ellipseI) {
    this->setLength(ellipseI.getL());
    this->setMatrixPhase(0);
    this->setInclusions_T(ellipseI.getMicroInclusions());
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::setInclusions(
    const vector<smallShape::SpheroPolyhedronInc<DIM>>& spheroPolyhedrons_, Point<DIM> L) {
    this->setLength(L);
    this->setMatrixPhase(0);
    setInnerShapes(spheroPolyhedrons_);
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::setInclusions(const Rectangle<DIM>& rect) {
    this->setLength(rect.getL());
    this->setMatrixPhase(0);
    vector < smallShape::ConvexPolyhedronInc<DIM>> thePolyhedrons = {
            smallShape::Rectangle<DIM>(1, create_array<DIM>(0.), rect.recL) };
    setInnerShapes(thePolyhedrons);
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::addLayer(
    const vector<Identifier>& identifiers,
    const vector<PhaseType>& newPhase, const vector<double>& width) {
    addLayer(smallShape::auxi_layerInstructions::buildInstructionVector(identifiers, newPhase, width));
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::addLayer(const vector<Identifier>& identifiers,
    PhaseType newPhase, double width) {
    vector<PhaseType> newPhaseList(identifiers.size(), newPhase);
    vector<double> widthList(identifiers.size(), width);
    addLayer(smallShape::auxi_layerInstructions::buildInstructionVector(identifiers, newPhaseList, widthList));
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::addLayer(
    vector<smallShape::LayerInstructions> layersToAdd) {
    checkAddLayer(layersToAdd);
    applyOnAllInclusions([&layersToAdd](auto& vectorInclusions) {
        auxi_MultiInclusions::addLayer_T(layersToAdd, vectorInclusions);
        });
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::setInclusions(
    const SpheroPolyhedronInclusions<DIM>& spheroPolyhedrons_) {
    setInclusions(spheroPolyhedrons_.getMicroInclusions(), spheroPolyhedrons_.getL());
}

template<unsigned short DIM>
inline bool MultiInclusions<DIM>::checkAddLayer(
    vector<smallShape::LayerInstructions> layersToAdd) const {
    ///
    auto errorMessage = [](Identifier identifier) {
        cerr << __PRETTY_FUNCTION__ << endl;
        cerr << "Identifier = " << identifier << endl;
        throw runtime_error("Unknown identifier!");
        };
    ///
    auto allIdentifiers = getAllIdentifiers();
    sort(allIdentifiers.begin(), allIdentifiers.end());
    sort(layersToAdd.begin(), layersToAdd.end());
    size_t j = 0;
    for (const auto& layer : layersToAdd) {
        while (true) {
            if (j > allIdentifiers.size()) {
                errorMessage(layer.identifier);
            }
            if (layer.identifier < allIdentifiers[j]) {
                errorMessage(layer.identifier);
            } else if (layer.identifier == allIdentifiers[j]) {
                break;
            } else {
                j++;
            }
        }
    }
    return true;
}

template<unsigned short DIM>
inline vector<Identifier> MultiInclusions<DIM>::getAllIdentifiers() const {
    vector<Identifier> allId{ };
    applyOnAllInclusions([&allId](const auto& vectorInclusions) {
        for (const auto& sph : vectorInclusions) {
            allId.push_back(sph.identifier);
        }
        });
    assert(allId == auxi_MultiInclusions::getAllIdentifiers(sphereInc.size() + polyhedrons.size() + ellipseInc.size() + spheroPolyhedrons.size()));
    return allId;
}

template<unsigned short DIM>
inline void MultiInclusions<DIM>::changePhase(const vector<Identifier>& identifiers,
    const vector<PhaseType>& newPhase) {
    vector <smallShape::LayerInstructions> layerInstructions = smallShape::auxi_layerInstructions::buildInstructionVector(identifiers, newPhase);
    applyOnAllInclusions([&layerInstructions](auto& vectorInclusions) {
        auxi_MultiInclusions::changePhase_T(layerInstructions, vectorInclusions);
        });
}

template<unsigned short DIM>
inline vector<Identifier> MultiInclusions<DIM>::getIdentifiers(
    vector<PhaseType> phases) const {
    vector<Identifier>identifiersAll{};
    applyOnAllInclusions([&identifiersAll, &phases](const auto& vectorInclusions) {
        auto identifiers1 = auxi_MultiInclusions::getIdentifiers(phases, vectorInclusions);
        std::copy(identifiers1.begin(), identifiers1.end(), std::back_inserter(identifiersAll));
        });
    return identifiersAll;
}

template<unsigned short DIM>
inline vector<Point<DIM> > MultiInclusions<DIM>::getAllCenters() const {
    vector < Point<DIM> > allId{ };
    applyOnAllInclusions([&allId](const auto& vectorInclusions) {
        for (const auto& sph : vectorInclusions) {
            allId.push_back(sph.center);
        }
        });
    return allId;
}

template<unsigned short DIM>
inline vector<PhaseType> MultiInclusions<DIM>::getAllPhases() const {
    std::set<PhaseType> allPhases{};
    applyOnAllInclusions([&allPhases](auto& vectorInclusions) {
        for (const auto& pol : vectorInclusions) {
            for (const auto j : pol.getLayerPhases()) {
                allPhases.insert(j);
            }
        }
        });
    if (allPhases.size() == 0 or sphereInc.size() != 0 or ellipseInc.size() != 0 or spheroPolyhedrons.size() != 0) {
        allPhases.insert(matrixPhase);
    }
    vector<PhaseType> result{ };
    copy(allPhases.begin(), allPhases.end(), back_inserter(result));
    return result;
}

template<unsigned short DIM>
inline vector<PhaseType> MultiInclusions<DIM>::getAllPhases(size_t layer_index) const {
    std::set<PhaseType> allPhases;
    applyOnAllInclusions([&allPhases, &layer_index](auto& vectorInclusions) {
        for (const auto& pol : vectorInclusions) {
            if (layer_index < pol.getNbOfLayers())
                allPhases.insert(pol.getPhaseGraphical(layer_index));
        }
        });
    //!
    if (layer_index > 0 and (allPhases.size() == 0 or sphereInc.size() != 0 or ellipseInc.size() != 0 or spheroPolyhedrons.size() != 0)) {
        allPhases.insert(matrixPhase);
    }
    vector<PhaseType> result{ };
    copy(allPhases.begin(), allPhases.end(), back_inserter(result));
    return result;
}

template<unsigned short DIM>
template<class LAMBDA_FUNCTION>
inline void merope::MultiInclusions<DIM>::applyOnAllInclusions(LAMBDA_FUNCTION f) {
    f(polyhedrons);
    f(sphereInc);
    f(ellipseInc);
    f(spheroPolyhedrons);
}

template<unsigned short DIM>
template<class LAMBDA_FUNCTION>
inline void merope::MultiInclusions<DIM>::applyOnAllInclusions(LAMBDA_FUNCTION f) const {
    f(polyhedrons);
    f(sphereInc);
    f(ellipseInc);
    f(spheroPolyhedrons);
}

//!
inline vector<Identifier> auxi_MultiInclusions::getAllIdentifiers(
    size_t NbOfSeeds) {
    vector < Identifier > result = { };
    for (size_t i = 0; i < NbOfSeeds; i++) {
        result.push_back(i);
    }
    return result;
}


template<class INCLUSIONVECTOR>
inline void auxi_MultiInclusions::changePhase_T(
    vector<smallShape::LayerInstructions> layerInstructions, INCLUSIONVECTOR& inclusions) {
    using C = typename INCLUSIONVECTOR::value_type;
    auto pointerInclusions = sortInclusionAndInstructions < C
    >(inclusions, layerInstructions);
    auto instruction = [](C* inclusion, smallShape::LayerInstructions layerInstruction) {
        inclusion->getPhaseGraphical(0) = layerInstruction.phase;
        };
    applyLayerInstruction_T < C
    >(layerInstructions, pointerInclusions, instruction);
}

template<class INCLUSIONVECTOR>
inline void auxi_MultiInclusions::addLayer_T(
    vector<smallShape::LayerInstructions> layerInstructions, INCLUSIONVECTOR& inclusions) {
    using C = typename INCLUSIONVECTOR::value_type;
    auto pointerInclusions = sortInclusionAndInstructions < C
    >(inclusions, layerInstructions);
    auto instruction = [](C* inclusion, smallShape::LayerInstructions layerInstruction) {
        inclusion->pushLayer(layerInstruction.width, layerInstruction.phase);
        };
    applyLayerInstruction_T<C>(layerInstructions, pointerInclusions, instruction);
}

template<class C>
inline void auxi_MultiInclusions::applyLayerInstruction_T(
    const vector<smallShape::LayerInstructions>& layersInstructions,
    vector<C*>& sortedPointerInclusionList,
    Instruction<C> applyInstruction) {
    size_t i_inc = 0;
    size_t i_lay = 0;
    while (i_inc < sortedPointerInclusionList.size()
        and i_lay < layersInstructions.size()) {
        const auto& instruction = layersInstructions[i_lay];
        auto& poly = sortedPointerInclusionList[i_inc];
        if (instruction.identifier == poly->identifier) {
            applyInstruction(poly, instruction);
            //poly->pushLayer(layerInstruction.width, layerInstruction.phase);
            i_inc++;
        } else if (instruction.identifier > poly->identifier) {
            i_inc++;
        } else {
            i_lay++;
        }
    }
}

template<class C>
inline vector<C*> auxi_MultiInclusions::sortInclusionAndInstructions(
    vector<C>& inclusionVector,
    vector<smallShape::LayerInstructions>& layerInstructions) {
    vector<C*> pointerInclusions{ }; // for sorting
    for (auto& inclusion : inclusionVector) {
        pointerInclusions.push_back(&inclusion);
    }
    vector<Identifier> index_id{ };
    sort(layerInstructions.begin(), layerInstructions.end(),
        [](const auto& l1, const auto& l2) {
            return l1.identifier < l2.identifier;
        });
    sort(pointerInclusions.begin(), pointerInclusions.end(),
        [](const auto& i1, const auto& i2) {
            return i1->identifier < i2->identifier;
        });
    return pointerInclusions;
}

template<class INCLUSIONVECTOR>
inline vector<Identifier> auxi_MultiInclusions::getIdentifiers(
    vector<PhaseType> phases, const INCLUSIONVECTOR& inclusionList) {
    //static_assert(std::is_base_of<smallShape::MicroInclusion<DIM>,C>::value);
    auto indicesInclusion = auxi_MultiInclusions::getAllIdentifiers(inclusionList.size());
    vector<Identifier> selectedIndices = {};
    //
    auto fun1 = [&inclusionList = std::as_const(inclusionList)](auto index) {
        return inclusionList[index].getPhaseGraphical(0);
        };
    auto fun2 = [](const PhaseType phase) {
        return phase;
        };
    auxi_function::extract_list(indicesInclusion.begin(), indicesInclusion.end(),
        phases.begin(), phases.end(),
        std::back_inserter(selectedIndices), fun1, fun2);
    //
    vector < Identifier > identifiers{ };
    std::transform(selectedIndices.begin(), selectedIndices.end(),
        std::back_inserter(identifiers), [&inclusionList = std::as_const(inclusionList)](auto i) {
            return inclusionList[i].identifier;
        });
    return identifiers;
}

} // namespace merope


#endif /* MULTIINCLUSIONS_IXX_ */
