//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_INCLUSIONANDVOXEL_IXX_
#define GRID_INCLUSIONANDVOXEL_IXX_

#include "../MeropeNamespace.hxx"


namespace merope {

template<class VOXEL_TYPE>
inline void vox::inclusionAndVoxel::fillVoxel(VOXEL_TYPE& voxelData, const VOXEL_TYPE& phases2Include) {
    static_assert(std::is_same<VOXEL_TYPE, VTK_PHASE>::value or std::is_same<VOXEL_TYPE, vox::VoxelPhaseFrac>::value);
    //
    if constexpr (std::is_same<VOXEL_TYPE, VTK_PHASE>::value) {
        voxelData = phases2Include;
    }
    else if constexpr (std::is_same<VOXEL_TYPE, vox::VoxelPhaseFrac>::value) {
        for (const auto& phfv : phases2Include) {
            voxelData.push_back(phfv);
        }
    }
}

template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
inline bool vox::inclusionAndVoxel::onlyOnePhaseInsideVoxel(const INCLUSION& microInclusion,
    const Point<DIM>& centerVoxel, const double& halfDiagVoxel,
    VOXEL_TYPE& phases2Include) {
    //
    for (size_t indexLayer = 0; indexLayer < microInclusion.getNbOfLayers(); indexLayer++) {
        if (microInclusion.guaranteeInside(centerVoxel - microInclusion.center, halfDiagVoxel, indexLayer)) {
            if (indexLayer == microInclusion.getNbOfLayers() - 1 or microInclusion.guaranteeOutside(centerVoxel - microInclusion.center, halfDiagVoxel, indexLayer)) {
                if constexpr (std::is_same<VOXEL_TYPE, vox::VoxelPhaseFrac>::value) {
                    vox::inclusionAndVoxel::fillVoxel<VOXEL_TYPE>(phases2Include, VoxelPhaseFrac{ SinglePhaseFrac(microInclusion.getPhaseForVoxellation(indexLayer), 1.) });
                }
                else if constexpr (std::is_same<VOXEL_TYPE, vox::VTK_PHASE>::value) {
                    vox::inclusionAndVoxel::fillVoxel<VOXEL_TYPE>(phases2Include, microInclusion.getPhaseForVoxellation(indexLayer));
                }
                return true;
            }
            else {
                return false;
            }
        }
    }
    return false;
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
    if constexpr (std::is_same<VOXEL_TYPE, vox::VTK_PHASE>::value) {
        Point<DIM> centerVoxel = vox::auxi::center <DIM>(indexVoxel, dx);
        auto layerIndex = microInclusion.whichLayer(centerVoxel);
        if (layerIndex >= 0) {
            phases2Include = microInclusion.getPhaseForVoxellation(layerIndex);
            return true;
        }
    }
    else if constexpr (std::is_same<VOXEL_TYPE, vox::VoxelPhaseFrac>::value) {
        Point<DIM> centerPoly_to_origVoxel = vox::auxi::origin <DIM>(indexVoxel, dx) - microInclusion.center;
        phases2Include = inclusionAndVoxel::computeAllFracVol<DIM, INCLUSION, VOXEL_TYPE>(microInclusion, centerPoly_to_origVoxel, dx, halfDiagVoxel);
        return true;
    }
    return false; // no intersection
}


template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
inline VOXEL_TYPE vox::inclusionAndVoxel::computeAllFracVol(
    const INCLUSION& inclusion, const Point<DIM>& centerPoly_to_origVoxel, const Point<DIM>& dx, const double& halfDiagVoxel) {
    static_assert(std::is_same<VOXEL_TYPE, vox::VoxelPhaseFrac>::value);
    static_assert(std::is_same<INCLUSION, smallShape::ConvexPolyhedronInc<DIM>>::value or std::is_same<INCLUSION, smallShape::SphereInc<DIM>>::value
        or std::is_same<INCLUSION, smallShape::EllipseInc<DIM>>::value or std::is_same<INCLUSION, smallShape::SpheroPolyhedronInc<DIM>>::value);

    VOXEL_TYPE phases2Include{ };
    if constexpr (std::is_same<INCLUSION, smallShape::ConvexPolyhedronInc<DIM>>::value) {
        const auto& innerInclusions = inclusion.getInnerInclusions();
        for (size_t layer_i = 0; layer_i < innerInclusions.size(); layer_i++) {
            double volFrac = geomTools::fracVolIntersection<DIM>(innerInclusions[layer_i].faces, centerPoly_to_origVoxel, dx);
            if (volFrac < geomTools::EPSILON) {
                break;
            }
            phases2Include.push_back(SinglePhaseFrac(inclusion.getPhaseForVoxellation(layer_i), volFrac));
        }
    }
    else if constexpr (std::is_same<INCLUSION, smallShape::SpheroPolyhedronInc<DIM>>::value) {
        const auto& innerInclusions = inclusion.getInnerInclusions();
        for (size_t layer_i = 0; layer_i < innerInclusions.size(); layer_i++) {
            double volFrac = geomTools::fracVolIntersection<DIM>(innerInclusions[layer_i], centerPoly_to_origVoxel, dx);
            if (volFrac < geomTools::EPSILON) {
                break;
            }
            phases2Include.push_back(SinglePhaseFrac(inclusion.getPhaseForVoxellation(layer_i), volFrac));
        }
    }
    else if constexpr (std::is_same<INCLUSION, smallShape::SphereInc<DIM>>::value) {
        if (not inclusion.guaranteeOutside(centerPoly_to_origVoxel, 2 * halfDiagVoxel)) { // otherwise, sure that the voxel does not intersect the sphere
            Point<DIM> normal = centerPoly_to_origVoxel + 0.5 * dx; // computed wrt to the center of the voxel
            geomTools::renormalize<DIM>(normal);
            HalfSpace<DIM> tangentPlane(normal, inclusion.getInnerInclusions()[0].radius);
            for (size_t layer_i = 0; layer_i < inclusion.getNbOfLayers();
                layer_i++) {
                tangentPlane.c() = inclusion.getInnerInclusions()[layer_i].radius;
                double volFrac = geomTools::fracVolIntersection<DIM>(tangentPlane, centerPoly_to_origVoxel, dx);
                if (volFrac < geomTools::EPSILON) {
                    break;
                }
                phases2Include.push_back(SinglePhaseFrac(inclusion.getPhaseForVoxellation(layer_i), volFrac));
            }
        }
    }
    else if constexpr (std::is_same<INCLUSION, smallShape::EllipseInc<DIM>>::value) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("To be programmed");
    }
    // dispatches the volume fractions between many phases
    for (size_t j = 0; j + 1 < phases2Include.size(); j++) {
        phases2Include[j].fracVol -= phases2Include[j + 1].fracVol;
    }
    return phases2Include;
}

} // namespace merope


#endif /* GRID_INCLUSIONANDVOXEL_IXX_ */
