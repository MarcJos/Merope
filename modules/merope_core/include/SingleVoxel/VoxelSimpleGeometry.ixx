//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../../Geometry/include/VolumeInCube.hxx"
#include "Nan_like.hxx"

namespace merope {
namespace vox {


template<unsigned short DIM, class BasicType>
double VoxelWithGeometry<DIM, BasicType>::remove_empty_parts() {
    // remove empty parts
    double volume = 0;
    vector<BasicType> new_phases{};
    vector<vector<HalfSpace<DIM>>> new_halfSpaces{};
    if (this->phases.size() > 0) {
        vector<bool> remove_poly(this->phases.size(), false);
        for (long i = 0; i < this->phases.size(); i++) {
            double vol = geomTools::volume_in_cube(this->halfSpaces[i]);
            if (vol > 0) {
                new_phases.push_back(this->phases[i]);
                new_halfSpaces.push_back(this->halfSpaces[i]);
                volume += vol;
            }
        }
        std::swap(new_phases, this->phases);
        std::swap(new_halfSpaces, this->halfSpaces);
    }
    return volume;
}

template<unsigned short DIM, class BasicType>
void VoxelWithGeometry<DIM, BasicType>::postProcess(
    const MatrixPhaseHolder<BasicType>& matrixPhaseHolder) {
    this->remove_empty_parts();
    if (matrixPhaseHolder.is_there_matrix()) {
        this->setMatrixPhase(matrixPhaseHolder.getMatrixPhase());
    }
}

template<unsigned short DIM, class BasicType>
void VoxelWithGeometry<DIM, BasicType>::setCombination2(const VoxelWithGeometry<DIM, BasicType>& voxel_0,
    const VoxelWithGeometry<DIM, BasicType>& voxel_1,
    std::function<BasicType(BasicType, BasicType)> func,
    std::function<bool(BasicType, BasicType)> func_filter) {
    auto intersections = compute_intersections(voxel_0, voxel_1);
    this->phases.resize(intersections.phases.size());
    this->halfSpaces.resize(intersections.halfSpaces.size());
    size_t j_index = 0;
    // transform the (phase_0, phase_1) of each inclusion into a func(phase_0, phase_1)
    auto updatePhase = [&](auto newPhase, size_t i) {
        this->phases[j_index] = newPhase;
        this->halfSpaces[j_index] = intersections.halfSpaces[i];
        j_index++;
        };
    auto insertPhase2 = [&](auto phi1, auto phi2, size_t i) {
        if (func_filter(phi1, phi2)) {
            auto newPhase = func(phi1, phi2);
            updatePhase(newPhase, i);
        }
        };
    auto insertPhase1 = [&](auto phi1, size_t i) {
        updatePhase(phi1, i);
        };
    //
    for (size_t i = 0; i < intersections.phases.size(); i++) {
        const auto& my_phases = intersections.phases[i];
        if (is_nan_like(get<0>(my_phases))) {
            if (voxel_0.is_there_matrix()) {
                insertPhase2(voxel_0.getMatrixPhase(), get<1>(my_phases), i);
            } else {
                insertPhase1(get<1>(my_phases), i);
            }
        } else if (is_nan_like(get<1>(my_phases))) {
            if (voxel_1.is_there_matrix()) {
                insertPhase2(get<0>(my_phases), voxel_1.getMatrixPhase(), i);
            } else {
                insertPhase1(get<0>(my_phases), i);
            }
        } else {
            insertPhase2(get<0>(my_phases), get<1>(my_phases), i);
        }
    }
    this->phases.resize(j_index);
    this->halfSpaces.resize(j_index);
    //
    // set matrix phase
    if (voxel_0.is_there_matrix() and voxel_1.is_there_matrix()) {
        if (func_filter(voxel_0.getMatrixPhase(), voxel_1.getMatrixPhase())) {
            this->setMatrixPhase(func(voxel_0.getMatrixPhase(), voxel_1.getMatrixPhase()));
        }
    }
    this->remove_empty_parts();
}


template<unsigned short DIM>
vector<vector<HalfSpace<DIM>>> compute_removal(const vector<HalfSpace<DIM>>& polyhedron,
    const vector<HalfSpace<DIM>>& polyhedron_to_remove) {
    vector<vector<HalfSpace<DIM>>> result{};
    if (polyhedron_to_remove.size() == 0) {
        return {};
    } else {
        for (size_t i = 0; i < polyhedron_to_remove.size(); i++) {
            auto poly = polyhedron;
            for (size_t j = 0; j < i + 1; j++) {
                auto hf = polyhedron_to_remove[j];
                if (j == i) {
                    hf.vec_force_definition() = RenormPoint<DIM>(-hf.vec());
                    hf.c() = -hf.c();
                }
                poly.push_back(hf);
            }
            result.push_back(poly);
        }
        return result;
    }
}

template<unsigned short DIM, class BasicType>
VoxelWithGeometry<DIM, std::pair<BasicType, BasicType>> compute_intersections(
    const VoxelWithGeometry<DIM, BasicType>& voxel_0,
    const VoxelWithGeometry<DIM, BasicType>& voxel_1) {
    VoxelWithGeometry<DIM, std::pair<BasicType, BasicType>> result{};
    auto add_phase = [&result](const auto& list_hf, auto phi_0, auto phi_1) {
        for (const auto& hf : list_hf) {
            result.halfSpaces.push_back(hf);
            result.phases.push_back(std::pair<BasicType, BasicType>(phi_0, phi_1));
        }
        };

    for (size_t i_0 = 0; i_0 < voxel_0.phases.size(); i_0++) {
        for (size_t i_1 = 0; i_1 < voxel_1.phases.size(); i_1++) {
            auto phi_0 = voxel_0.phases[i_0];
            auto phi_1 = voxel_1.phases[i_1];
            auto hf_0 = voxel_0.halfSpaces[i_0];
            auto hf_1 = voxel_1.halfSpaces[i_1];
            ////////////////////
            vector<vector<HalfSpace<DIM>>> list_hf_0 = compute_removal(hf_0, hf_1);
            add_phase(list_hf_0, phi_0, make_nan_like<BasicType>());
            ////////////////////
            vector<vector<HalfSpace<DIM>>> list_hf_1 = compute_removal(hf_1, hf_0);
            add_phase(list_hf_1, make_nan_like<BasicType>(), phi_1);
            ////////////////////
            auto intersection = hf_0;
            for (size_t j = 0; j < hf_1.size(); j++) {
                intersection.push_back(hf_1[j]);
            }
            result.halfSpaces.push_back(intersection);
            result.phases.push_back(std::pair<BasicType, BasicType>({ phi_0, phi_1 }));
        }
    }
    // special case no polyhedra inside 1 volume
    if (voxel_1.phases.size() == 0) {
        result.halfSpaces = voxel_0.halfSpaces;
        for (const auto& phi_0 : voxel_0.phases) {
            result.phases.push_back(std::pair<BasicType, BasicType>({ phi_0,  make_nan_like<BasicType>() }));
        }
    }
    // special case no polyhedra inside 1 volume
    if (voxel_0.phases.size() == 0) {
        result.halfSpaces = voxel_1.halfSpaces;
        for (const auto& phi_1 : voxel_1.phases) {
            result.phases.push_back(std::pair<BasicType, BasicType>({ make_nan_like<BasicType>(), phi_1 }));
        }
    }
    return result;
}

template<unsigned short DIM, class BasicType>
void VoxelWithGeometry<DIM, BasicType>::merge_with(const VoxelWithGeometry<DIM, BasicType>& anotherVoxel,
    const std::function<BasicType(BasicType, BasicType)>& ruleIntersection) {
    if (anotherVoxel.is_empty()) {
        // nothing
    } else if (this->is_empty()) {
        this->halfSpaces = anotherVoxel.halfSpaces;
        this->phases = anotherVoxel.phases;
    } else {
        auto result_intersections = compute_intersections(*this, anotherVoxel);
        this->halfSpaces = result_intersections.halfSpaces;
        this->phases.resize(result_intersections.phases.size());
        for (size_t i = 0; i < result_intersections.phases.size(); i++) {
            const auto& my_phases = result_intersections.phases[i];
            if (is_nan_like(get<0>(my_phases))) {
                this->phases[i] = get<1>(my_phases);
            } else if (is_nan_like(get<1>(my_phases))) {
                this->phases[i] = get<0>(my_phases);
            } else {
                this->phases[i] = ruleIntersection(get<0>(my_phases), get<1>(my_phases));
            }
        }
    }
}

template<unsigned short DIM, class BasicType>
void VoxelWithGeometry<DIM, BasicType>::merge_with(const VoxelWithGeometry<DIM, BasicType>& anotherVoxel) {
    // matrix phase
    if (anotherVoxel.is_there_matrix()) {
        if (this->is_there_matrix()) {
            Merope_assert(anotherVoxel.getMatrixPhase() == anotherVoxel.getMatrixPhase(), "merging two different voxels only possible if same matrix phase");
        }
        this->setMatrixPhase(anotherVoxel.getMatrixPhase());
    }
    //
    if (anotherVoxel.is_empty()) {
        // nothing
    } else if (this->is_empty()) {
        this->halfSpaces = anotherVoxel.halfSpaces;
        this->phases = anotherVoxel.phases;
    } else {
        for (size_t i = 0; i < anotherVoxel.phases.size(); i++) {
            this->phases.push_back(anotherVoxel.phases[i]);
            this->halfSpaces.push_back(anotherVoxel.halfSpaces[i]);
        }
    }
}

}  // namespace  vox
}  // namespace  merope