//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../../../GenericMerope/StdHeaders.hxx"
#include "../../../GenericMerope/MeropeNamespace.hxx"
#include "../../../Geometry/include/AmbiantSpace.hxx"

#include "RecurStructure.hxx"
#include "MatrixPhaseHolder.hxx"

namespace merope {
namespace vox {

template<unsigned short DIM, class BasicType>
class VoxelWithGeometry : public MatrixPhaseHolder<BasicType> {
    //! @brief voxel [0, 1]^DIM with inside phases delimitated by HalfSpaces
public:
    //! constructors
    VoxelWithGeometry() : MatrixPhaseHolder<BasicType>(), phases{}, halfSpaces{} {}
    VoxelWithGeometry(BasicType phase) : MatrixPhaseHolder<BasicType>(), phases{ phase }, halfSpaces{ {} } {}

    //! @brief insert the data of another voxel. Assume to make a unique of inclusions
    //! @param anotherVoxel 
    //! @param ruleIntersection
    void merge_with(const VoxelWithGeometry<DIM, BasicType>& anotherVoxel,
        const std::function<BasicType(BasicType, BasicType)>& ruleIntersection);
    //! assume no intersection
    void merge_with(const VoxelWithGeometry<DIM, BasicType>& anotherVoxel);

    //! @brief : set a pure phase inside the voxel
    void set_pure(BasicType purePhase) {
        this->phases = { purePhase };
        this->halfSpaces = { {} };
    }

    //! @return whether the voxel is empty (=no inner polyhedron inside)
    bool is_empty() const { return phases.size() == 0; }

    //! @brief remove the empty components of voxels
    double remove_empty_parts();
    //! @brief remove the empty components of voxels
    //! + insert matrix phase 
    void postProcess(const MatrixPhaseHolder<BasicType>& matrixPhaseHolder);

    //! @brief : make a combination of 2 voxels though a function
    void setCombination2(const VoxelWithGeometry<DIM, BasicType>& voxel_0,
        const VoxelWithGeometry<DIM, BasicType>& voxel_1,
        std::function<BasicType(BasicType, BasicType)> func,
        std::function<bool(BasicType, BasicType)> func_filter = [](BasicType, BasicType) {return true;});

    //! @brief 
    using BASIC_TYPE = BasicType;
public:
    vector<BasicType> phases;
    vector<vector<HalfSpace<DIM>>> halfSpaces;
};


//! @return the list of non-intersecting polyhedra equal to polyhedron \ polyhedron_to_remove
//! @param polyhedron 
//! @param polyhedron_to_remove 
template<unsigned short DIM>
vector<vector<HalfSpace<DIM>>> compute_removal(const vector<HalfSpace<DIM>>& polyhedron,
    const vector<HalfSpace<DIM>>& polyhedron_to_remove);

//! @return compute the resulting intersection of 2 different geometries
//! @warning : assume each geometry contains non-intersecting polyhedra
//! @param voxel_0 
//! @param voxel_1 
template<unsigned short DIM, class BasicType>
VoxelWithGeometry<DIM, std::pair<BasicType, BasicType>> compute_intersections(
    const VoxelWithGeometry<DIM, BasicType>& voxel_0,
    const VoxelWithGeometry<DIM, BasicType>& voxel_1);

}  // namespace  vox
}  // namespace  merope

#include "VoxelSimpleGeometry.ixx"
