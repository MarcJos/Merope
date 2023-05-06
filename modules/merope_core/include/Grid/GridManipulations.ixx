//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_GRIDMANIPULATIONS_IXX_
#define GRID_GRIDMANIPULATIONS_IXX_


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM, class VOXEL_TYPE, class FUNCTION>
CartesianGrid<DIM, VOXEL_TYPE> gridAuxi::combineGridFunction(const CartesianGrid<DIM, VOXEL_TYPE>& grid1,
    const CartesianGrid<DIM, VOXEL_TYPE>& grid2, const FUNCTION  func) {
    CartesianGrid<DIM, VOXEL_TYPE> result = grid1;
#pragma omp parallel for firstprivate(func)
    for (size_t i = 0; i < grid1.size(); i++) {
        result[i] = gridAuxi::combineVoxelFunc(grid1[i], grid2[i], func);
    }
    return result;
}

template<unsigned short DIM, class VOXEL_TYPE>
CartesianGrid<DIM, VOXEL_TYPE> gridAuxi::combineGridMask(const CartesianGrid<DIM, VOXEL_TYPE>& grid1,
    const CartesianGrid<DIM, VOXEL_TYPE>& grid2,
    const CartesianGrid<DIM, VOXEL_TYPE>& mask) {
    CartesianGrid<DIM, VOXEL_TYPE> result = grid1;
#pragma omp parallel for
    for (size_t i = 0; i < grid1.size(); i++) {
        result[i] = gridAuxi::combineVoxelMask(grid1[i], grid2[i], mask[i]);
    }
    return result;
}

template<class VOXEL_TYPE>
VOXEL_TYPE gridAuxi::combineVoxelMask(const VOXEL_TYPE& phf1, const VOXEL_TYPE& phf2, const VOXEL_TYPE& mask) {
    static_assert(std::is_same<VOXEL_TYPE, vox::OutputFormat<VoxelRule::Average>>::value
        or std::is_arithmetic<VOXEL_TYPE>::value);
    //
    if constexpr (std::is_same<VOXEL_TYPE, vox::OutputFormat<VoxelRule::Average>>::value) { // case of PhaseFracVol
        array<double, 2> proportions = vox::gridAuxi::translateMask(mask);
        if (proportions[1] < geomTools::EPSILON) {
            return phf1;
        }
        else if (proportions[0] < geomTools::EPSILON) {
            return phf2;
        }
        else {
            VoxelPhaseFrac result{};
            for (auto phf : phf1) {
                phf.fracVol *= proportions[0];
                result.push_back(phf);
            }
            for (auto phf : phf2) {
                phf.fracVol *= proportions[1];
                result.push_back(phf);
            }
            return result;
        }
    }
    else if constexpr (std::is_arithmetic<VOXEL_TYPE>::value) { // case of numeric type
        if (mask > 0) {
            return phf2;
        }
        else {
            return phf1;
        }
    }
}

template<class VOXEL_TYPE, class FUNCTION>
VOXEL_TYPE gridAuxi::combineVoxelFunc(const VOXEL_TYPE& voxphf1, const VOXEL_TYPE& voxphf2, const FUNCTION& func) {
    static_assert(std::is_same<VOXEL_TYPE, vox::OutputFormat<VoxelRule::Average>>::value
        or std::is_arithmetic<VOXEL_TYPE>::value);
    //
    if constexpr (std::is_same<VOXEL_TYPE, vox::OutputFormat<VoxelRule::Average>>::value) { // case of PhaseFracVol
        VoxelPhaseFrac result{};
        for (const auto& pf1 : voxphf1) {
            for (const auto& pf2 : voxphf2) {
                result.push_back(
                    SinglePhaseFrac(func(pf1.phase, pf2.phase), pf1.fracVol * pf2.fracVol));
            }
        }
        result.merge();
        return result;
    }
    else if constexpr (std::is_arithmetic<VOXEL_TYPE>::value) { // case of numeric type
        return func(voxphf1, voxphf2);
    }
}

} // namespace vox
} // namespace merope

#endif /* GRID_GRIDMANIPULATIONS_IXX_ */
