//! Copyright : see license.txt
//!
//! \brief


#pragma once

#include "../../AlgoPacking/include/SphereContainer.hxx"
#include "../../Voxellation/VoxRecurStructure.hxx"


#include "../../VTKinout/VTK_adapter.hxx"


namespace merope {
namespace vox {

// Voxellation_data
template<unsigned short DIM>
inline Voxellation_data<DIM>::Voxellation_data(array<double, DIM> L) :
    GridParameters<DIM>(create_array<DIM, size_t>(1), L),
    voxelRule{ vox::VoxelRule::Center },
    homogRule{ homogenization::Rule::Smallest },
    structure{ nullptr }, fieldStructure{ nullptr },
    pureCoeffs{}, pureCoeffsIsSet{ false }{
}

template<unsigned short DIM>
inline Voxellation_data<DIM>::Voxellation_data(const Structure<DIM>& structure_) :
    Voxellation_data(structure_.getL()) {
    structure.reset(new Structure<DIM>(structure_));
}

template<unsigned short DIM>
inline Voxellation_data<DIM>::Voxellation_data(const FieldStructure<DIM>& fieldStructure_) :
    Voxellation_data(fieldStructure_.getL()) {
    fieldStructure.reset(new FieldStructure<DIM>(fieldStructure_));
}

template<unsigned short DIM>
void Voxellation_data<DIM>::setDiscretization(array<size_t, DIM> nbNodes_, array<size_t, DIM> nMin_, array<size_t, DIM> nMax_) {
    this->set_nbNodes_L(nbNodes_, this->getL());
    this->setSubGridIndices(nMin_, nMax_);
}

template<unsigned short DIM>
template<VoxelRule VOXEL_RULE>
void Voxellation_data<DIM>::check_voxel_rule(const char* pretty_function) const {
    if (VOXEL_RULE != this->voxelRule) {
        cerr << pretty_function << endl;
        throw invalid_argument("Incorrect VoxelRule");
    }
}

template<unsigned short DIM>
template<homogenization::Rule HOMOG_RULE>
void Voxellation_data<DIM>::check_homog_rule(const char* pretty_function) const {
    if (HOMOG_RULE != this->homogRule) {
        cerr << pretty_function << endl;
        throw invalid_argument("Incorrect HomogenizationRule");
    }
}

template<unsigned short DIM>
template<class STRUCTURE>
void Voxellation_data<DIM>::check_structure(const char* pretty_function) const {
    if ((is_same_v<STRUCTURE, Structure<DIM>> and not structure) or (is_same_v<STRUCTURE, FieldStructure<DIM>> and not fieldStructure)) {
        cerr << pretty_function << endl;
        throw invalid_argument("No structure is loaded");
    } else {
        if (not structure and not fieldStructure) {
            cerr << pretty_function << endl;
            throw invalid_argument("No structure is loaded");
        }
    }
}

template<unsigned short DIM>
void Voxellation_data<DIM>::setPureCoeffs(vector<double> pureCoeffs_) {
    pureCoeffsIsSet = true;
    pureCoeffs = pureCoeffs_;
}

template<unsigned short DIM>
void Voxellation_data<DIM>::checkPureCoeffs() {
    check_structure<Structure<DIM>>(__PRETTY_FUNCTION__);
    if (not pureCoeffsIsSet) {
        pureCoeffsIsSet = true;
        vector<PhaseType> phases{};
        phases = structure->getAllPhases();
        PhaseType maxPhase = *(max_element(phases.begin(), phases.end()));
        if (*(min_element(phases.begin(), phases.end())) < 0) {
            throw runtime_error("One phase is negative!");
        }
        vector<double> defaultPureCoeffs(maxPhase + 1);
        for (long i = 0; i < maxPhase + 1; i++) {
            defaultPureCoeffs[i] = i;
        }
        setPureCoeffs(defaultPureCoeffs);
    }
    // error verification
    if (this->homogRule == homogenization::Rule::Reuss) {
        for (auto pc : this->pureCoeffs) {
            Merope_assert(pc > 0,
                "A negative value is forbidden!");
        }
    }
}

// Voxellation

template<unsigned short DIM>
inline GridField<DIM> Voxellation<DIM>::getField() {
    this->template check_structure<>(__PRETTY_FUNCTION__);
    if (this->fieldStructure) {
        return computeField(*(this->fieldStructure));
    } else if (this->structure) {
        return computeField(*(this->structure));
    } else {
        Merope_assert(false, "Unexpected");
    }
}

template<unsigned short DIM>
template<class STRUCTURE>
inline GridField<DIM> Voxellation<DIM>::computeField(const STRUCTURE& structure_) {
    if (this->voxelRule == VoxelRule::Center) {
        if constexpr (is_same_v<STRUCTURE, FieldStructure<DIM>>) {
            return voxellizer::transformStructIntoGrid<DIM>(structure_, GridParameters<DIM>(*this), vox::VoxelPolicy<VoxelRule::Center, true>());
        } else {
            Merope_error_not_done();
        }
    } else {
        auto phaseValueField = voxellizer::transformStructIntoGrid<DIM>(structure_, GridParameters<DIM>(*this), vox::VoxelPolicy<VoxelRule::Average, true>());
        return applyHomogRule(phaseValueField);
    }
}

template<unsigned short DIM>
inline void Voxellation<DIM>::printFieldFile(string fileVTK) {
    auto field = this->getField();
    merope::vtk_adapter::printVTK<double>(field, fileVTK, "Value");
}

template<unsigned short DIM>
void Voxellation<DIM>::printFile(string fileVTK, string fileCoeff) {
    //
    this->template check_structure<>(__PRETTY_FUNCTION__);
    if (this->structure) this->checkPureCoeffs();
    //
    vector<double> coefficients{};
    auto gridPhase = computePurePhaseGrid(coefficients);
    // really print
    vtk_adapter::printVTK<unsigned short>(gridPhase, fileVTK);
    vtk_adapter::printCoeffs<DIM>(coefficients, fileCoeff);
}

template<unsigned short DIM>
string Voxellation<DIM>::print() {
    this->template check_structure<>(__PRETTY_FUNCTION__);
    if (this->structure) this->checkPureCoeffs();
    //
    vector<double> coefficients{};
    auto gridPhase = computePurePhaseGrid(coefficients);
    Grid_VER grid_ver{};
    grid_ver.set_L<DIM>(this->getL());
    grid_ver.set_Nb<DIM>(this->getNbNodeBigGrid());
    grid_ver.fromGridPhase<DIM>(gridPhase, coefficients);
    grid_ver.print();
    return "Printed in std output";
}

template<unsigned short DIM>
template<VoxelRule VOXEL_RULE>
CartesianGrid<DIM, composite::OutputFormat<VOXEL_RULE, DIM, PhaseType>> Voxellation<DIM>::computeCartesianGrid() {
    this-> template check_structure<Structure<DIM>>(__PRETTY_FUNCTION__);
    this-> template check_voxel_rule<VOXEL_RULE>(__PRETTY_FUNCTION__);
    return  voxellizer::transformStructIntoGrid<DIM>(*(this->structure), GridParameters<DIM>(*this), vox::VoxelPolicy<VOXEL_RULE, true>());
}

template<unsigned short DIM>
inline vector<composite::stl_format_Iso<PhaseType>> Voxellation<DIM>::computePhaseGrid() {
    unique_ptr<CartesianGrid<DIM, composite::Iso<PhaseType>>> res;
    if (this->voxelRule == VoxelRule::Average) {
        res.reset(new CartesianGrid<DIM, composite::Iso<PhaseType>>(computeCartesianGrid<VoxelRule::Average>()));
    } else {
        Merope_assert(false, "Incoherent voxel rule");
    }
    return convertGrid::linearize(convertGrid::convert_to_stl_format(*res));
}

template<unsigned short DIM>
vector<composite::stl_format_AnIso<DIM, PhaseType>> Voxellation<DIM>::computeCompositeGrid() {
    return convertGrid::linearize(convertGrid::convert_to_stl_format(
        computeCartesianGrid<VoxelRule::Laminate>()));
}

template<unsigned short DIM>
CartesianGrid<DIM, PhaseType> Voxellation<DIM>::computePurePhaseGrid(vector<double>& coefficients_) {
    coefficients_ = {};
    if (this->structure and this->voxelRule == VoxelRule::Center) {
        coefficients_ = this->pureCoeffs;
        auto gridPhase = voxellizer::transformStructIntoGrid<DIM>(*(this->structure), GridParameters<DIM>(*this), vox::VoxelPolicy<VoxelRule::Center, true>());
        convertGrid::renormalizeWithCoefficients(gridPhase, coefficients_);
        return gridPhase;
    } else {
        auto gridPhase = convertGrid::fromFieldToPhase(this->getField(), coefficients_);
        convertGrid::renormalizeWithCoefficients(gridPhase, coefficients_);
        return gridPhase;
    }
}

template<unsigned short DIM>
template<class VOXEL_TYPE>
GridField<DIM> Voxellation<DIM>::applyHomogRule(const CartesianGrid<DIM, VOXEL_TYPE>& phaseFracVol) {
    if (std::is_same_v<VOXEL_TYPE, composite::Iso<PhaseType>>) {
        this->checkPureCoeffs();
    }
    /////////////////
    if (this->homogRule == homogenization::Rule::Voigt) {
        return applyHomogRule_T<homogenization::Rule::Voigt>(phaseFracVol);
    } else if (this->homogRule == homogenization::Rule::Reuss) {
        return applyHomogRule_T<homogenization::Rule::Reuss>(phaseFracVol);
    } else if (this->homogRule == homogenization::Rule::Largest) {
        return applyHomogRule_T<homogenization::Rule::Largest>(phaseFracVol);
    } else if (this->homogRule == homogenization::Rule::Smallest) {
        return applyHomogRule_T<homogenization::Rule::Smallest>(phaseFracVol);
    } else {
        Merope_assert(false, "invalid_argument HomogRule");
    }
}

template<unsigned short DIM>
template<homogenization::Rule HOMOG_RULE, class VOXEL_TYPE>
vox::GridField<DIM> Voxellation<DIM>::applyHomogRule_T(const vox::CartesianGrid<DIM, VOXEL_TYPE>& phaseFracVol) {
    if constexpr (is_same_v<VOXEL_TYPE, composite::Iso<double>>) {
        return vox::voxellizer::applyHomogRule_T<HOMOG_RULE, DIM>(phaseFracVol);
    } else if constexpr (is_same_v<VOXEL_TYPE, composite::Iso<PhaseType>>) {
        return vox::voxellizer::applyHomogRule_T<HOMOG_RULE, DIM>(phaseFracVol, this->pureCoeffs);
    } else {
        Merope_assert(false, "Incorrect VOXEL_TYPE");
    }
}

}  // namespace vox
}  // namespace merope




