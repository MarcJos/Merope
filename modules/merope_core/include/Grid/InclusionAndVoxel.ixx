//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../MeropeNamespace.hxx"
#include "../Geometry/GeomTools.hxx"


namespace merope {

template<class VOXEL_TYPE>
inline void vox::inclusionAndVoxel::fillVoxel(VOXEL_TYPE& voxelData, const VOXEL_TYPE& phases2Include) {
    static_assert(std::is_same_v<VOXEL_TYPE, PhaseType>
        or vox::composite::is_Iso<VOXEL_TYPE>
        or vox::composite::is_AnIso<VOXEL_TYPE>);
    //
    if constexpr (vox::composite::is_Pure<VOXEL_TYPE>) {
        voxelData = phases2Include;
    } else if constexpr (vox::composite::is_Iso<VOXEL_TYPE> or vox::composite::is_AnIso<VOXEL_TYPE>) {
        for (const auto& phfv : phases2Include) {
            voxelData.push_back(phfv);
        }
    }
}

template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
bool vox::inclusionAndVoxel::fillVoxel(
    const INCLUSION& microInclusion,
    const DiscPoint<DIM>& indexVoxel, const Point<DIM>& dx,
    const double& halfDiagVoxel, VOXEL_TYPE& voxelData) {
    VOXEL_TYPE phases2Include;
    bool isThereIntersection = vox::inclusionAndVoxel::phasesInsideVoxel<DIM, INCLUSION, VOXEL_TYPE>(microInclusion, indexVoxel, dx, halfDiagVoxel, phases2Include);
    if (isThereIntersection) {
        vox::inclusionAndVoxel::fillVoxel<VOXEL_TYPE>(voxelData, phases2Include);
    }
    return isThereIntersection;
}

template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
inline bool vox::inclusionAndVoxel::phasesInsideVoxel(
    const INCLUSION& microInclusion, const DiscPoint<DIM>& indexVoxel, const Point<DIM>& dx, const double& halfDiagVoxel, VOXEL_TYPE& phases2Include) {
    if constexpr (std::is_same_v<VOXEL_TYPE, vox::PhaseType>) {
        Point<DIM> centerVoxel = vox::auxi::center <DIM>(indexVoxel, dx);
        auto layerIndex = microInclusion.whichLayer(centerVoxel);
        if (layerIndex >= 0) {
            phases2Include = microInclusion.getPhaseForVoxellation(layerIndex);
            return true;
        }
    } else if constexpr (composite::is_Iso<VOXEL_TYPE> or composite::is_AnIso<VOXEL_TYPE>) {
        Point<DIM> centerPoly_to_origVoxel = vox::auxi::origin <DIM>(indexVoxel, dx) - microInclusion.center;
        phases2Include = inclusionAndVoxel::computeAllFracVol<DIM, INCLUSION, VOXEL_TYPE>(microInclusion, centerPoly_to_origVoxel, dx, halfDiagVoxel);
        return true;
    }
    return false;  // no intersection
}


template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
inline VOXEL_TYPE vox::inclusionAndVoxel::computeAllFracVol(
    const INCLUSION& inclusion, const Point<DIM>& centerPoly_to_origVoxel, const Point<DIM>& dx, const double& halfDiagVoxel) {
    static_assert(composite::is_Iso<VOXEL_TYPE> or composite::is_AnIso<VOXEL_TYPE>);
    static_assert(smallShape::IsInc<DIM, INCLUSION>);

    VOXEL_TYPE phases2Include{ };
    Point<DIM> vector_intersection{};

    auto insert_volfrac = [&](size_t layer_i, double volFrac) {
        if (volFrac < geomTools::EPSILON) {
            return false;
        } else {
            if constexpr (composite::is_Iso<VOXEL_TYPE>) {
                phases2Include.push_back(typename VOXEL_TYPE::PHASE_FRAC(inclusion.getPhaseForVoxellation(layer_i), volFrac));
            } else if constexpr (composite::is_AnIso<VOXEL_TYPE>) {
                phases2Include.push_back(typename VOXEL_TYPE::PHASE_FRAC(inclusion.getPhaseForVoxellation(layer_i), volFrac, -vector_intersection));
            }
            return true;
        }
        };

    if constexpr (
        std::is_same_v<INCLUSION, smallShape::ConvexPolyhedronInc<DIM>>
        or std::is_same_v<INCLUSION, smallShape::SpheroPolyhedronInc<DIM>>
        or (DIM == 3 and std::is_same_v<INCLUSION, smallShape::CylinderInc<3>>)) {
        const auto& innerInclusions = inclusion.getInnerInclusions();
        for (size_t layer_i = 0; layer_i < innerInclusions.size(); layer_i++) {
            double volFrac = geomTools::fracVolIntersection<DIM>(innerInclusions[layer_i], centerPoly_to_origVoxel, dx, vector_intersection);
            if (not insert_volfrac(layer_i, volFrac)) {
                break;
            }
        }
    } else if constexpr (std::is_same_v<INCLUSION, smallShape::SphereInc<DIM>>) {
        if (not inclusion.guaranteeOutside(centerPoly_to_origVoxel, 2 * halfDiagVoxel)) {  // otherwise, sure that the voxel does not intersect the sphere
            Point<DIM> normal = centerPoly_to_origVoxel + 0.5 * dx;  // computed wrt to the center of the voxel
            geomTools::renormalize<DIM>(normal);
            HalfSpace<DIM> tangentPlane(normal, inclusion.getInnerInclusions()[0].radius);
            for (size_t layer_i = 0; layer_i < inclusion.getNbOfLayers();
                layer_i++) {
                tangentPlane.c() = inclusion.getInnerInclusions()[layer_i].radius;
                vector_intersection = tangentPlane.vec();
                double volFrac = geomTools::fracVolIntersection<DIM>(tangentPlane, centerPoly_to_origVoxel, dx);
                if (not insert_volfrac(layer_i, volFrac)) {
                    break;
                }
            }
        }
    } else if constexpr (std::is_same_v<INCLUSION, smallShape::EllipseInc<DIM>>) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("To be programmed");
    }
    // dispatches the volume fractions between many phases
    for (size_t j = 0; j + 1 < phases2Include.size(); j++) {
        phases2Include[j].fracVol -= phases2Include[j + 1].fracVol;
    }
    return phases2Include;
}

}  // namespace merope



