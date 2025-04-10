//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../Geometry/include/GeomTools.hxx"
#include "../../../Geometry/include/GeomTools_1.hxx"


namespace merope {

template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE, class VOXEL_POLICY>
bool vox::inclusionAndVoxel::fillVoxel(
    const INCLUSION& microInclusion,
    const DiscPoint<DIM>& indexVoxel, const Point<DIM>& dx,
    const double& halfDiagVoxel, VOXEL_TYPE& voxelData,
    const VOXEL_POLICY& voxelPolicy) {
    VOXEL_TYPE phases2Include;
    bool isThereIntersection = vox::inclusionAndVoxel::phasesInsideVoxel<DIM, INCLUSION, VOXEL_TYPE>(microInclusion, indexVoxel, dx, halfDiagVoxel, phases2Include);
    if (isThereIntersection) {
        vox::composite::voxel_filler<DIM>(voxelData, phases2Include, voxelPolicy);
    }
    return isThereIntersection;
}

template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
inline bool vox::inclusionAndVoxel::phasesInsideVoxel(
    const INCLUSION& microInclusion, const DiscPoint<DIM>& indexVoxel, const Point<DIM>& dx, const double& halfDiagVoxel, VOXEL_TYPE& phases2Include) {
    if constexpr (std::is_same_v<VOXEL_TYPE, PhaseType>) {
        Point<DIM> centerVoxel = vox::auxi::center <DIM>(indexVoxel, dx);
        auto layerIndex = microInclusion.whichLayer(centerVoxel);
        if (layerIndex >= 0) {
            phases2Include = microInclusion.getPhaseForVoxellation(layerIndex);
            return true;
        }
    } else if constexpr (composite::is_Iso<VOXEL_TYPE> or composite::is_AnIso<VOXEL_TYPE>) {
        Point<DIM> centerPoly_to_origVoxel = vox::auxi::origin <DIM>(indexVoxel, dx) - microInclusion.center;
        phases2Include = inclusionAndVoxel::computeAllFracVol<DIM, INCLUSION, VOXEL_TYPE>(microInclusion, centerPoly_to_origVoxel, dx, halfDiagVoxel);
        return (phases2Include.size() > 0);
    } else if constexpr (composite::is_PolyGeom<VOXEL_TYPE>) {
        Point<DIM> centerPoly_to_origVoxel = vox::auxi::origin <DIM>(indexVoxel, dx) - microInclusion.center;
        phases2Include = inclusionAndVoxel::computeSimpleGeometry<DIM, INCLUSION, VOXEL_TYPE>(microInclusion, centerPoly_to_origVoxel, dx);
        return (not phases2Include.is_empty());
    } else {
        Merope_static_error(VOXEL_TYPE, "Incorrect voxel Rule");
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
        Merope_error_not_done();
    }
    // dispatches the volume fractions between many phases
    for (size_t j = 0; j + 1 < phases2Include.size(); j++) {
        phases2Include[j].fracVol -= phases2Include[j + 1].fracVol;
    }
    return phases2Include;
}

template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
VOXEL_TYPE vox::inclusionAndVoxel::computeSimpleGeometry(const INCLUSION& inclusion, const Point<DIM>& centerPoly_to_origVoxel, const Point<DIM>& dx) {
    static_assert(composite::is_PolyGeom<VOXEL_TYPE>);
    static_assert(smallShape::IsInc<DIM, INCLUSION>);
    VoxelWithGeometry<DIM, typename VOXEL_TYPE::BASIC_TYPE> res{};
    //
    vector<vector<HalfSpace<DIM>>> list_of_tangent_planes(inclusion.getNbOfLayers());
    size_t max_layer = inclusion.getNbOfLayers();
    for (size_t layer_i = 0; layer_i < inclusion.getNbOfLayers(); layer_i++) {
        bool is_empty = false;
        std::tie(is_empty, list_of_tangent_planes[layer_i])
            = get_list_tangentPlanes<DIM, INCLUSION>(inclusion, centerPoly_to_origVoxel, dx, layer_i);
        if (is_empty) {
            max_layer = layer_i;
            break;
        }
    }
    //
    for (size_t layer_i = 0; layer_i < max_layer; layer_i++) {
        auto phase = inclusion.getPhaseForVoxellation(layer_i);
        if (layer_i + 1 < max_layer) {
            auto union_of_poly_layer = vox::compute_removal<DIM>(list_of_tangent_planes[layer_i], list_of_tangent_planes[layer_i + 1]);
            for (size_t j = 0; j < union_of_poly_layer.size(); j++) {
                res.phases.push_back(phase);
                res.halfSpaces.push_back(union_of_poly_layer[j]);
            }
        } else { // layer_i + 1 == max_layer
            res.phases.push_back(phase);
            res.halfSpaces.push_back(list_of_tangent_planes[layer_i]);
        }
    }
    //
    VOXEL_TYPE finalResult(std::move(res));
    return finalResult;
}

template<unsigned short DIM, class INCLUSION>
tuple<bool, vector<HalfSpace<DIM>>> vox::inclusionAndVoxel::get_list_tangentPlanes(const INCLUSION& inclusion, const Point<DIM>& centerPoly_to_origVoxel, const Point<DIM>& dx, size_t layer_i) {
    bool is_empty = false;
    //
    Point<DIM> inverse_dx = dx;
    for (auto& coord : inverse_dx) {
        coord = 1. / coord;
    }
    //
    vector<HalfSpace<DIM>> list_tangentPlanes = {};
    auto renormalize_and_add_tangentPlane = [&centerPoly_to_origVoxel, &list_tangentPlanes, &is_empty, &inverse_dx](HalfSpace<DIM> tangentPlane) {
        // renormalization
        tangentPlane.c() -= geomTools::prodScal<DIM>(tangentPlane.vec(), centerPoly_to_origVoxel);
        linearTransform::proceed<DIM>(tangentPlane, inverse_dx);
        // 
        double vol_frac = geomTools::fracVolIntersection<DIM>(tangentPlane);
        if (vol_frac < geomTools::EPSILON) {
            is_empty = true;
            // the intersection is empty
        } else if (vol_frac > 1 - geomTools::EPSILON) {
            // this plane is irrelevant, for it does not intersect the voxel
        } else {
            list_tangentPlanes.push_back(tangentPlane);
        }
        };
    //
    const auto& innerInclusions_i = inclusion.getInnerInclusions()[layer_i];
    auto hs_intersecting = geomTools::make_HalfSpace_Intersecting_Voxel<DIM>(innerInclusions_i,
        centerPoly_to_origVoxel, dx);

    while (not hs_intersecting.is_empty()) {
        renormalize_and_add_tangentPlane(*(hs_intersecting.get()));
        hs_intersecting.next();
    }
    return { is_empty, list_tangentPlanes };
}

}  // namespace merope



