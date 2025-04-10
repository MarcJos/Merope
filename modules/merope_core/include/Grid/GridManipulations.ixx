//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


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
    if constexpr (composite::is_Iso<VOXEL_TYPE>) {  // case of PhaseFracVol
        array<double, 2> proportions = vox::gridAuxi::translateMask(mask);
        if (proportions[1] < geomTools::EPSILON) {
            return phf1;
        } else if (proportions[0] < geomTools::EPSILON) {
            return phf2;
        } else {
            VOXEL_TYPE result{};
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
    } else if constexpr (composite::is_Pure<VOXEL_TYPE>) {  // case of numeric type
        if (mask > 0) {
            return phf2;
        } else {
            return phf1;
        }
    } else if constexpr (composite::is_AnIso<VOXEL_TYPE>) {
        Merope_error_not_done();
    } else if constexpr (composite::is_PolyGeom<VOXEL_TYPE>) {
        VOXEL_TYPE vox1, vox2;
        vox1.setCombination2(phf1, mask,
            [](auto phase, auto /*phaseMask*/) {return phase;},
            [](auto /*phase*/, auto phaseMask) {return phaseMask <= 0;});
        vox2.setCombination2(phf2, mask,
            [](auto phase, auto /*phaseMask*/) {return phase;},
            [](auto /*phase*/, auto phaseMask) {return phaseMask > 0;});
        vox2.merge_with(vox1);
        return vox2;
    } else {
        Merope_static_error(VOXEL_TYPE, "Unknown voxel rule");
    }
}

template<class VOXEL_TYPE, class FUNCTION>
VOXEL_TYPE gridAuxi::combineVoxelFunc(const VOXEL_TYPE& voxphf1, const VOXEL_TYPE& voxphf2, const FUNCTION& func) {
    if constexpr (composite::is_Pure<VOXEL_TYPE>) {  // case of numeric type
        return func(voxphf1, voxphf2);
    } else if constexpr (composite::is_Iso<VOXEL_TYPE>
        or composite::is_AnIso<VOXEL_TYPE>) {  // case of PhaseFracVol
        VOXEL_TYPE result{};
        for (const auto& pf1 : voxphf1) {
            for (const auto& pf2 : voxphf2) {
                auto pf3 = pf1;
                pf3.phase = func(pf1.phase, pf2.phase);
                pf3.fracVol = pf1.fracVol * pf2.fracVol;
                if constexpr (composite::is_AnIso<VOXEL_TYPE>) {
                    if (voxphf1.size() == 1) { pf3.normal = pf2.normal; }
                }
                result.push_back(pf3);
            }
        }
        result.merge();
        return result;
    } else if constexpr (composite::is_PolyGeom<VOXEL_TYPE>) {
        VOXEL_TYPE result{};
        result.setCombination2(voxphf1, voxphf2, func);
        return result;
    } else {
        Merope_static_error(VOXEL_TYPE, "Incorrect voxel type");
    }
}

}  // namespace vox
}  // namespace merope


