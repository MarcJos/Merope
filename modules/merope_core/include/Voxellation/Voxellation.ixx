//! Copyright : see license.txt
//!
//! \brief


#ifndef VOXELLATION_IXX_
#define VOXELLATION_IXX_

#include "../../../AlgoPacking/src/SphereContainer.hxx"
#include "../Voxellation/VoxRecurStructure.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

// Voxellation_data
template<unsigned short DIM>
inline Voxellation_data<DIM>::Voxellation_data(array<double, DIM> L):
    PreVoxellation<DIM>(L),
    homogRule{ homogenization::Rule::Smallest },
    structure{ nullptr }, fieldStructure{ nullptr },
    coefficients{}, pureCoeffs{}, pureCoeffsIsSet{ false }{
}

template<unsigned short DIM>
inline Voxellation_data<DIM>::Voxellation_data(const Structure<DIM>& structure_):
    Voxellation_data(structure_.getL()) {
    structure.reset(new Structure<DIM>(structure_));
}

template<unsigned short DIM>
inline Voxellation_data<DIM>::Voxellation_data(const FieldStructure<DIM>& fieldStructure_):
    Voxellation_data(fieldStructure_.getL()) {
    fieldStructure.reset(new FieldStructure<DIM>(fieldStructure_));
}

template<unsigned short DIM>
void Voxellation_data<DIM>::setDiscretization(array<size_t, DIM> nbNodes_, array<size_t, DIM> nMin_, array<size_t, DIM> nMax_) {
    this->setGridNL(nbNodes_, this->getL());
    this->setSubGridIndices(nMin_, nMax_);
}

template<unsigned short DIM>
bool Voxellation_data<DIM>::checkStructure(const char* pretty_function) {
    if (not this->structure and not this->fieldStructure) {
        cerr << pretty_function << endl;
        throw invalid_argument("No structure is loaded");
    }
    return true;
}

template<unsigned short DIM>
void vox::Voxellation_data<DIM>::setPureCoeffs(vector<double> pureCoeffs_) {
    pureCoeffsIsSet = true;
    pureCoeffs = pureCoeffs_;
    coefficients = pureCoeffs_;
}

template<unsigned short DIM>
void Voxellation_data<DIM>::checkPureCoeffs() {
    assert(this->structure);
    if (not pureCoeffsIsSet) {
        pureCoeffsIsSet = true;
        vector<PhaseType> phases{};
        phases = structure->getAllPhases();
        PhaseType maxPhase = *(max_element(phases.begin(), phases.end()));
        if (*(min_element(phases.begin(), phases.end())) < 0) {
            throw runtime_error("One phase is negative!");
        }
        vector<double> defaultPureCoeffs(maxPhase + 1);
        for (size_t i = 0; i < maxPhase + 1; i++) {
            defaultPureCoeffs[i] = i;
        }
        setPureCoeffs(defaultPureCoeffs);
    }
    // error verification
    if (this->homogRule == homogenization::Rule::Reuss) {
        for (auto pc : this->pureCoeffs) {
            if (pc <= 0) {
                cerr << __PRETTY_FUNCTION__ << endl;
                throw runtime_error("A negative value is forbidden!");
            }
        }
    }
}

template<unsigned short DIM>
void Voxellation_data<DIM>::printFile(string fileVTK, string fileCoeff) const {
    VTKstream vtkstream(fileVTK.c_str());
    this->grid.toVTKCELL(vtkstream);
    printCoeffs(fileCoeff);
}

template<unsigned short DIM>
void Voxellation_data<DIM>::printCoeffs(string fileCoeff) const {
    ofstream ost(fileCoeff);
    for (auto c : coefficients) {
        ost << c << endl;
    }
}

// Voxellation

template<unsigned short DIM>
void Voxellation<DIM>::proceed(array<size_t, DIM> nbNodes_) {
    proceed(nbNodes_, create_array<DIM, size_t>(0), nbNodes_);
}

template<unsigned short DIM>
void Voxellation<DIM>::proceed(array<size_t, DIM> nbNodes_, array<size_t, DIM> nMin_, array<size_t, DIM> nMax_) {
    this->setDiscretization(nbNodes_, nMin_, nMax_);
    //
    this->checkStructure(__PRETTY_FUNCTION__);
    if (this->structure) this->checkPureCoeffs();
    //
    auto gridPhase = getPhaseGrid();
    this->grid.template fromGridPhase<DIM>(gridPhase, this->coefficients);
}

template<unsigned short DIM>
inline GridField<DIM> Voxellation<DIM>::getField() {
    this->checkStructure(__PRETTY_FUNCTION__);
    if (this->fieldStructure) {
        return getField_FieldStructure();
    }
    else if (this->structure) {
        return getField_Structure();
    }
    else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected");
    }
}

template<unsigned short DIM>
inline GridField<DIM> Voxellation<DIM>::getField_Structure() {
    assert(this->structure);
    if (this->voxelRule == vox::VoxelRule::Center) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("The Center Voxel Rule cannot give rise to a field");
    }
    auto phaseFracVol = getPhaseFracField<VoxelRule::Average>();
    return applyHomogRule(phaseFracVol);
}

template<unsigned short DIM>
inline GridField<DIM> Voxellation<DIM>::getField_FieldStructure() {
    if (not this->fieldStructure) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("PhaseFracField only makes sense for fieldStructure");
    }
    //
    return VoxStructure<DIM, double, CartesianField<DIM>, double>(*(this->fieldStructure), *this)();
}

template<unsigned short DIM>
template<VoxelRule VOXEL_RULE>
inline CartesianGrid<DIM, OutputFormat<VOXEL_RULE>> Voxellation<DIM>::getPhaseFracField() {
    if (not this->structure) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("PhaseFracField only makes sense for structures");
    }
    //
    return VoxStructure<DIM, vox::OutputFormat<VOXEL_RULE>, MultiInclusions<DIM>, PhaseType>(*(this->structure), *this)();
}

template<unsigned short DIM>
inline vector<vector<tuple<vox::VTK_PHASE, double>>> vox::Voxellation<DIM>::computePhaseGrid(array<size_t, DIM> nbNodes_) {
    auto nMin_ = create_array<DIM, size_t>(0), nMax_ = nbNodes_;
    this->setDiscretization(nbNodes_, nMin_, nMax_);
    ////////
    if (not this->structure) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw invalid_argument("No Structure");
    }
    ////////
    auto VOXEL_RULE = this->voxelRule;
    if (VOXEL_RULE == vox::VoxelRule::Average) {
        return convertGrid::fromCartesianToVector(getPhaseFracField<vox::VoxelRule::Average>());
    }
    else if (VOXEL_RULE == vox::VoxelRule::Center) {
        return convertGrid::fromCartesianToVector(convertGrid::fromPhaseToFracVol(
            getPhaseFracField<vox::VoxelRule::Center>()));
    }
    else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected");
    }
}

template<unsigned short DIM>
vox::CartesianGrid<DIM, vox::VTK_PHASE> Voxellation<DIM>::getPhaseGrid() {
    if (this->structure and this->voxelRule == vox::VoxelRule::Center) {
        return getPhaseFracField<vox::VoxelRule::Center>();
    }
    else {
        return vox::convertGrid::fromFieldToPhase(this->getField(), this->coefficients);
    }
}

template<unsigned short DIM>
GridField<DIM> Voxellation<DIM>::applyHomogRule(const CartesianGrid<DIM, VoxelPhaseFrac>& phaseFracVol) {
    this->checkPureCoeffs();
    /////////////////
    if (this->homogRule == homogenization::Rule::Voigt) {
        return applyHomogRule_T<homogenization::Rule::Voigt>(phaseFracVol);
    }
    else if (this->homogRule == homogenization::Rule::Reuss) {
        return applyHomogRule_T<homogenization::Rule::Reuss>(phaseFracVol);
    }
    else if (this->homogRule == homogenization::Rule::Largest) {
        return applyHomogRule_T<homogenization::Rule::Largest>(phaseFracVol);
    }
    else if (this->homogRule == homogenization::Rule::Smallest) {
        return applyHomogRule_T<homogenization::Rule::Smallest>(phaseFracVol);
    }
    else {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw invalid_argument("HomogRule");
    }
}

template<unsigned short DIM>
template<homogenization::Rule HOMOG_RULE>
GridField<DIM> Voxellation<DIM>::applyHomogRule_T(const CartesianGrid<DIM, VoxelPhaseFrac>& phaseFracVol) {
    auto copyPureCoeffs = this->pureCoeffs;
    CartesianGrid<DIM, double> gridField = vox::convertGrid::localConvert<DIM, double, vox::VoxelPhaseFrac>(phaseFracVol,
        [&copyPureCoeffs](const auto& phaseFracVolLoc) {
            auto input = gridAuxi::getTabCoeff(phaseFracVolLoc, copyPureCoeffs);
            return homogenization::homog<HOMOG_RULE>(input[0], input[1]);
        }
    );
    return gridField;
}

} // namespace vox
} // namespace merope



#endif /* VOXELLATION_IXX_ */
