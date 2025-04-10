//! Copyright : see license.txt
//!
//! \brief
//


#include "MicroInclusion/MicroInclusion.hxx"


namespace merope {

size_t smallShape::auxi_MicroInclusions::getIndexPhaseGraphical(size_t phaseIndex, size_t nbOfLayers) {
    Merope_assert((phaseIndex <= nbOfLayers),
        "Too large index.");
    if (phaseIndex == 0) {
        return nbOfLayers - 1;
    } else {
        return phaseIndex - 1;
    }
}

vector<smallShape::LayerInstructions> smallShape::auxi_layerInstructions::buildInstructionVector(
    vector<Identifier> identifier, vector<PhaseType> phase,
    vector<double> width) {
    if (identifier.size() == 0) {
        return vector<LayerInstructions> {};
    } else if (width.size() == 0) {
        return buildInstructionVector(identifier, phase, vector<double>(identifier.size(), 0.));
    }
    Merope_assert(identifier.size() == phase.size() and identifier.size() == width.size(),
        "Incompatible sizes");
    vector < LayerInstructions > layerInstructions{ };
    for (size_t i = 0; i < identifier.size(); i++) {
        layerInstructions.push_back(LayerInstructions{ identifier[i], phase[i],
                width[i] });
    }
    return layerInstructions;
}

}  // namespace merope

