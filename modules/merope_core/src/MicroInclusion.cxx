//! Copyright : see license.txt
//!
//! \brief
//


#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "MicroInclusion/MicroInclusion.hxx"


#include "MeropeNamespace.hxx"


namespace merope {

size_t smallShape::auxi_MicroInclusions::getIndexPhaseGraphical(size_t phaseIndex, size_t nbOfLayers) {
    if (phaseIndex > nbOfLayers) {
        cerr << __PRETTY_FUNCTION__;
        throw runtime_error("Too large index.");
    }
    if (phaseIndex == 0) {
        return nbOfLayers - 1;
    }
    else {
        return phaseIndex - 1;
    }
}

vector<smallShape::LayerInstructions> smallShape::auxi_layerInstructions::buildInstructionVector(
    vector<Identifier> identifier, vector<PhaseType> phase,
    vector<double> width) {
    if (identifier.size() == 0) {
        return vector<LayerInstructions> {};
    }
    else if (width.size() == 0) {
        return buildInstructionVector(identifier, phase, vector<double>(identifier.size(), 0.));
    }
    if (identifier.size() != phase.size() or identifier.size() != width.size()) {
        cerr << __PRETTY_FUNCTION__;
        throw runtime_error("Incompatible sizes");
    }
    vector < LayerInstructions > layerInstructions{ };
    for (size_t i = 0; i < identifier.size(); i++) {
        layerInstructions.push_back(LayerInstructions{ identifier[i], phase[i],
                width[i] });
    }
    return layerInstructions;
}

} // namespace merope

