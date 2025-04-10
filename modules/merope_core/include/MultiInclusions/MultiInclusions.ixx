//! Copyright : see license.txt
//!
//! \brief
//!

#pragma once


#include "../../../GenericTools/CPP_Functions.hxx"


namespace merope {
template<unsigned short DIM>
template<class C>
void MultiInclusions<DIM>::setInnerShapes(const C& inclusions) {
    // reinitialize the shapes
    apply_on_all([](auto& vectorInclusions) {
        vectorInclusions = {};
        });
    // fill in the correct shape
    using INCLUSION_TYPE = typename C::value_type;
    vector<INCLUSION_TYPE>& vectorShapes = this->template get<INCLUSION_TYPE>();
    vectorShapes = inclusions;
    //!
    stable_sort(vectorShapes.begin(), vectorShapes.end(),
        [](const auto& poly1, const auto& poly2) {
            return poly1.identifier < poly2.identifier;
        });
    //! Identifier shall be the position of the multiInclusion
    for (size_t i = 0; i < inclusions.size(); i++) {
        vectorShapes[i].identifier = i;
    }
    laguerreTess = nullptr;
}

template<unsigned short DIM>
MultiInclusions<DIM>::MultiInclusions() :
    InsideTorus<DIM>(), SOA_type(), MatrixPhaseHolder<PhaseType>{}, laguerreTess{ nullptr }{
}

template<unsigned short DIM>
void MultiInclusions<DIM>::setInclusions(LaguerreTess<DIM> polyX) {
    this->setLength(polyX.getL());
    polyX.computeTessels();
    setInnerShapes(polyX.getMicroInclusions());
    this->laguerreTess.reset(new LaguerreTess<DIM>(polyX));
}

template<unsigned short DIM>
void MultiInclusions<DIM>::setInclusions(const PolyInclusions<DIM>& polyInc) {
    this->setLength(polyInc.getL());
    this->setMatrixPhase_if_not_present(0);
    setInnerShapes(polyInc.getMicroInclusions());
}

template<unsigned short DIM>
void MultiInclusions<DIM>::setInclusions(
    const SphereInclusions<DIM>& sphereI) {
    this->setLength(sphereI.getL());
    this->setMatrixPhase_if_not_present(0);

    vector<smallShape::SphereInc<DIM>> theSphereInc{};
    transform(sphereI.getSpheres().begin(), sphereI.getSpheres().end(), std::back_inserter(theSphereInc), [](const auto& sph) {
        return smallShape::SphereInc<DIM>(sph);
        });
    setInnerShapes(theSphereInc);
}

template<unsigned short DIM>
template<class ObjectInc>
void MultiInclusions<DIM>::setInclusions(const ObjectInclusions<DIM, ObjectInc>& objectInclusion) {
    setInclusions(objectInclusion.getMicroInclusions(), objectInclusion.getL());
}

template<unsigned short DIM>
template<class ObjectInc>
void MultiInclusions<DIM>::setInclusions(
    const vector<ObjectInc>& vectInclusions, Point<DIM> L) {
    this->setLength(L);
    this->setMatrixPhase_if_not_present(0);
    setInnerShapes(vectInclusions);
}

template<unsigned short DIM>
void MultiInclusions<DIM>::setInclusions(const Rectangle<DIM>& rect) {
    this->setLength(rect.getL());
    this->setMatrixPhase_if_not_present(0);
    vector<smallShape::ConvexPolyhedronInc<DIM>> thePolyhedrons = {
            smallShape::Rectangle<DIM>(1, create_array<DIM>(0.), rect.recL) };
    setInnerShapes(thePolyhedrons);
}

template<unsigned short DIM>
void MultiInclusions<DIM>::addLayer(
    const vector<Identifier>& identifiers,
    const vector<PhaseType>& newPhase, const vector<double>& width) {
    addLayer(smallShape::auxi_layerInstructions::buildInstructionVector(identifiers, newPhase, width));
}

template<unsigned short DIM>
void MultiInclusions<DIM>::addLayer(const vector<Identifier>& identifiers,
    PhaseType newPhase, double width) {
    vector<PhaseType> newPhaseList(identifiers.size(), newPhase);
    vector<double> widthList(identifiers.size(), width);
    addLayer(smallShape::auxi_layerInstructions::buildInstructionVector(identifiers, newPhaseList, widthList));
}

template<unsigned short DIM>
void MultiInclusions<DIM>::addLayer(
    vector<smallShape::LayerInstructions> layersToAdd) {
    checkAddLayer(layersToAdd);
    apply_on_all([&layersToAdd](auto& vectorInclusions) {
        auxi_MultiInclusions::addLayer_T(layersToAdd, vectorInclusions);
        });
}

template<unsigned short DIM>
void MultiInclusions<DIM>::enlarge(const vector<Identifier>& identifiers, const vector<double>& width) {
    Merope_warning(not isLaguerreTess(), "Laguerre tessellation : impossible to enlarge (no space between tessels)");
    auto instructions = smallShape::auxi_layerInstructions::buildInstructionVector(identifiers, vector<PhaseType>(width.size()), width);
    apply_on_all([&instructions](auto& vectorInclusions) {
        auxi_MultiInclusions::enlarge_T(instructions, vectorInclusions);
        });
}

template<unsigned short DIM>
bool MultiInclusions<DIM>::checkAddLayer(
    vector<smallShape::LayerInstructions> layersToAdd) const {
    ///
    auto errorMessage = [](Identifier identifier) {
        cerr << __PRETTY_FUNCTION__ << endl;
        cerr << "Identifier = " << identifier << endl;
        throw runtime_error("Unknown identifier!");
        };
    auto errorMessage_width = [](double width) {
        cerr << __PRETTY_FUNCTION__ << endl;
        cerr << "width = " << width << endl;
        throw runtime_error("Only positive width are allowed!");
        };
    ///
    auto allIdentifiers = getAllIdentifiers();
    sort(layersToAdd.begin(), layersToAdd.end());
    size_t j = 0;
    for (const auto& layer : layersToAdd) {
        if (layer.width < 0) {
            errorMessage_width(layer.width);
        }
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
vector<Identifier> MultiInclusions<DIM>::getAllIdentifiers() const {
    return auxi_MultiInclusions::getAllIdentifiers(this->get_total_size());
}

template<unsigned short DIM>
void MultiInclusions<DIM>::changePhase(const vector<Identifier>& identifiers,
    const vector<PhaseType>& newPhase) {
    vector <smallShape::LayerInstructions> layerInstructions = smallShape::auxi_layerInstructions::buildInstructionVector(identifiers, newPhase);
    apply_on_all([&layerInstructions](auto& vectorInclusions) {
        auxi_MultiInclusions::changePhase_T(layerInstructions, vectorInclusions);
        });
}

template<unsigned short DIM>
vector<Identifier> MultiInclusions<DIM>::getIdentifiers(
    vector<PhaseType> phases) const {
    vector<Identifier>identifiersAll{};
    apply_on_all([&identifiersAll, &phases](const auto& vectorInclusions) {
        auto identifiers1 = auxi_MultiInclusions::getIdentifiers(phases, vectorInclusions);
        std::copy(identifiers1.begin(), identifiers1.end(), std::back_inserter(identifiersAll));
        });
    return identifiersAll;
}

template<unsigned short DIM>
vector<Point<DIM> > MultiInclusions<DIM>::getAllCenters() const {
    vector < Point<DIM> > allId{ };
    apply_on_all([&allId](const auto& vectorInclusions) {
        for (const auto& sph : vectorInclusions) {
            allId.push_back(sph.center);
        }
        });
    return allId;
}

template<unsigned short DIM>
vector<PhaseType> MultiInclusions<DIM>::getAllPhases(size_t layer_index) const {
    auto test_and_fill_function = [&layer_index](auto& allPhases, const auto& pol) {
        if (layer_index < pol.getNbOfLayers())
            allPhases.insert(pol.getPhaseGraphical(layer_index));
        };
    const auto& this_to_ref = *this;
    auto insert_matrix_phase = [&layer_index, &this_to_ref](const auto&) {
        return false;
        };
    return getAllPhases(test_and_fill_function, insert_matrix_phase);
}

template<unsigned short DIM>
vector<PhaseType> MultiInclusions<DIM>::getAllPhases() const {
    auto test_and_fill_function = [](auto& allPhases, const auto& pol) {
        for (const auto j : pol.getLayerPhases()) {
            allPhases.insert(j);
        }
        };
    bool matrix_presence_ = this->is_there_matrix();
    auto insert_matrix_phase = [matrix_presence_](const auto& allPhases) {
        return allPhases.size() == 0 and matrix_presence_;
        };
    return getAllPhases(test_and_fill_function, insert_matrix_phase);
}

template<unsigned short DIM>
template<class TEST_AND_FILL_FUNCTION, class TEST_FUNCTION>
vector<PhaseType> MultiInclusions<DIM>::getAllPhases(TEST_AND_FILL_FUNCTION test_and_fill_function,
    TEST_FUNCTION insert_matrix_phase) const {
    std::set<PhaseType> allPhases;
    this->apply_on_all([&allPhases, &test_and_fill_function](auto& vectorInclusions) {
        for (const auto& pol : vectorInclusions) {
            test_and_fill_function(allPhases, pol);
        }
        });
    //!
    if (insert_matrix_phase(allPhases)) {
        allPhases.insert(getMatrixPhase());
    }
    vector<PhaseType> result{ };
    copy(allPhases.begin(), allPhases.end(), back_inserter(result));
    return result;
}

template<unsigned short DIM>
const LaguerreTess<DIM>& MultiInclusions<DIM>::getLaguerreTess() const {
    if (laguerreTess) {
        return *laguerreTess;
    } else {
        std::cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("No unerlying Laguerre tessellation!");
    }
}

//!
vector<Identifier> auxi_MultiInclusions::getAllIdentifiers(
    size_t NbOfSeeds) {
    vector < Identifier > result = { };
    for (Identifier i = 0; i < NbOfSeeds; i++) {
        result.push_back(i);
    }
    return result;
}

template<class INCLUSIONVECTOR>
void auxi_MultiInclusions::enlarge_T(vector<smallShape::LayerInstructions> all_instructions,
    INCLUSIONVECTOR& inclusions) {
    using C = typename INCLUSIONVECTOR::value_type;
    auto pointerInclusions = sortInclusionAndInstructions<C>(inclusions, all_instructions);
    auto instruction = [](C* inclusion, const smallShape::LayerInstructions& instruction_) {
        inclusion->enlarge(instruction_.width);
        };
    applyLayerInstruction_T<C>(all_instructions, pointerInclusions, instruction);
}

template<class INCLUSIONVECTOR>
void auxi_MultiInclusions::changePhase_T(
    vector<smallShape::LayerInstructions> layerInstructions, INCLUSIONVECTOR& inclusions) {
    using C = typename INCLUSIONVECTOR::value_type;
    auto pointerInclusions = sortInclusionAndInstructions<C>(inclusions, layerInstructions);
    auto instruction = [](C* inclusion, smallShape::LayerInstructions layerInstruction) {
        inclusion->getPhaseGraphical(0) = layerInstruction.phase;
        };
    applyLayerInstruction_T<C>(layerInstructions, pointerInclusions, instruction);
}

template<class INCLUSIONVECTOR>
void auxi_MultiInclusions::addLayer_T(
    vector<smallShape::LayerInstructions> layerInstructions, INCLUSIONVECTOR& inclusions) {
    using C = typename INCLUSIONVECTOR::value_type;
    auto pointerInclusions = sortInclusionAndInstructions<C>(inclusions, layerInstructions);
    auto instruction = [](C* inclusion, smallShape::LayerInstructions layerInstruction) {
        inclusion->pushLayer(layerInstruction.width, layerInstruction.phase);
        };
    applyLayerInstruction_T<C>(layerInstructions, pointerInclusions, instruction);
}

template<class C>
void auxi_MultiInclusions::applyLayerInstruction_T(
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
    vector<C*> pointerInclusions{ };  // for sorting
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
vector<Identifier> auxi_MultiInclusions::getIdentifiers(
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
    vector<Identifier> identifiers{ };
    std::transform(selectedIndices.begin(), selectedIndices.end(),
        std::back_inserter(identifiers), [&inclusionList = std::as_const(inclusionList)](auto i) {
            return inclusionList[i].identifier;
        });
    return identifiers;
}

}  // namespace merope



