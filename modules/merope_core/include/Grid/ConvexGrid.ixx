//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef GRID_CONVEXGRID_IXX_
#define GRID_CONVEXGRID_IXX_

#include "../MeropeNamespace.hxx"



namespace merope {
namespace vox {

// convexGrid<DIM>

template<unsigned short DIM>
inline auxi::ConvexGrid<DIM>::ConvexGrid(const array<array<long, 2>, DIM>& gridLimits_,
    const array<double, DIM>& dx_):
    dx{ dx_ }, halfDiagonal{ 0 },
    gridLimits{ gridLimits_ }, gridSize{},
    allFirstIndices{ getAllIndices_multi<DIM, DIM - 1, long>(gridLimits_) },
    sliceInstructions{} {
    for (size_t i = 0;i < DIM; i++) {
        gridSize[i] = gridLimits[i][1] - gridLimits[i][0];
    }
    sliceInstructions.resize(getNbFirstIndices());
    halfDiagonal = 0.5 * sqrt(geomTools::normeCarre<DIM>(dx));
    //
    assert(allFirstIndices.size() == sliceInstructions.size());
}



template<unsigned short DIM>
inline geomTools::Intersection_LineConvex vox::auxi::auxi_convexGrid::getLimitEllipse(
    array<double, 2>& dbLocalLimits, const Ellipse<DIM>& sphere,
    const Point<DIM - 1>& x1x2, double distance) {
    //
    cerr << __PRETTY_FUNCTION__ << endl;
    throw runtime_error("Should be programmed");
}


template<unsigned short DIM>
inline geomTools::Intersection_LineConvex vox::auxi::auxi_convexGrid::getLimitSphere(
    array<double, 2>& newLimits, const Sphere<DIM>& sphere, const Point<DIM - 1>& x1x2, double distance) {
    double heightWrtCenter = auxi_function::puissance<2>(sphere.radius + distance) - geomTools::normeCarre<DIM - 1>(x1x2);
    // first case : empty intersection
    if (heightWrtCenter < 0) {
        return geomTools::Intersection_LineConvex::Empty;
    }
    // second case : a segment
    heightWrtCenter = sqrt(heightWrtCenter);
    newLimits[0] = sphere.center[DIM - 1] - heightWrtCenter;
    newLimits[1] = sphere.center[DIM - 1] + heightWrtCenter;
    return geomTools::Intersection_LineConvex::Segment;
}

template<unsigned short DIM>
inline geomTools::Intersection_LineConvex vox::auxi::auxi_convexGrid::getLimitPolyhedron(
    array<double, 2>& dbLocalLimits, const ConvexPolyhedron<DIM>& polyhedron,
    const Point<DIM - 1>& x1x2, double distance) {
    double newLimit;
    dbLocalLimits = { -(0.5 * numeric_limits<double>::max()), 0.5 * numeric_limits<double>::max() };
    for (const auto& halfSpace : polyhedron.faces) {
        geomTools::Intersection_LineConvex typeIntersection = geomTools::computeIntersection(halfSpace, x1x2, newLimit, distance); // beware, 0 = center of polyhedron
        //
        if (typeIntersection == geomTools::Intersection_LineConvex::Minus) {
            dbLocalLimits[1] = min(dbLocalLimits[1], newLimit);
        }
        else if (typeIntersection == geomTools::Intersection_LineConvex::Plus) {
            dbLocalLimits[0] = max(dbLocalLimits[0], newLimit);
        }
        else if (typeIntersection == geomTools::Intersection_LineConvex::Empty) {
            return geomTools::Intersection_LineConvex::Empty;
        }
        // to avoid further tests
        if (dbLocalLimits[0] >= dbLocalLimits[1]) {
            return geomTools::Intersection_LineConvex::Empty;
        }
    }
    //
    dbLocalLimits[0] += polyhedron.center[DIM - 1];
    dbLocalLimits[1] += polyhedron.center[DIM - 1];
    return geomTools::Intersection_LineConvex::Segment;
}

template<unsigned short DIM>
inline geomTools::Intersection_LineConvex merope::vox::auxi::auxi_convexGrid::getLimitSpheroPolyhedron(
    array<double, 2>& newLimits,
    const SpheroPolyhedron<DIM>& spheroPolyhedron,
    const Point<DIM - 1>& x1x2, double distance) {
    if (distance <= 0) {
        return getLimitPolyhedron(newLimits, ConvexPolyhedron<DIM>(spheroPolyhedron.getCenter(), spheroPolyhedron.getInnerPolyhedron()), x1x2, distance);
    }
    else {
        return getLimitPolyhedron(newLimits, ConvexPolyhedron<DIM>(spheroPolyhedron.getCenter(), spheroPolyhedron.getOuterPolyhedron()), x1x2, distance);
    }
}


template<unsigned short DIM, class C>
inline bool vox::auxi::auxi_convexGrid::getLimits(const C& smallShape, const array<double, DIM>& center, const DiscPoint<DIM>& ij,
    const array<double, DIM>& dx, array<long, 2>& localLimits, const array<array<long, 2>, DIM>& gridLimits, double distance, VoxelType voxelType) {
    static_assert(std::is_same<C, sac_de_billes::ConvexPolyhedron<DIM>>::value or std::is_same<C, Sphere<DIM>>::value
        or std::is_same<C, Ellipse<DIM>>::value or std::is_same<C, SpheroPolyhedron<DIM>>::value);
    //
    double drift = -0.5 * dx[DIM - 1];
    array<double, DIM - 1> x1x2{};
    for (size_t i = 0; i < DIM - 1; i++) {
        x1x2[i] = (ij[i] + 0.5) * dx[i] - center[i];
    }
    //
    array<double, 2> dbLocalLimits;
    if constexpr (std::is_same<C, ConvexPolyhedron<DIM>>::value) {
        if (auxi_convexGrid::getLimitPolyhedron<DIM>(dbLocalLimits, smallShape, x1x2, distance) != geomTools::Intersection_LineConvex::Segment) {
            return false;
        }
    }
    else if constexpr (std::is_same<C, Sphere<DIM>>::value) {
        if (auxi_convexGrid::getLimitSphere<DIM>(dbLocalLimits, smallShape, x1x2, distance) != geomTools::Intersection_LineConvex::Segment) {
            return false;
        }
    }
    else if constexpr (std::is_same<C, Ellipse<DIM>>::value) {
        if (auxi_convexGrid::getLimitEllipse<DIM>(dbLocalLimits, smallShape, x1x2, distance) != geomTools::Intersection_LineConvex::Segment) {
            return false;
        }
    }
    else if constexpr (std::is_same<C, SpheroPolyhedron<DIM>>::value) {
        if (auxi_convexGrid::getLimitSpheroPolyhedron<DIM>(dbLocalLimits, smallShape, x1x2, distance) != geomTools::Intersection_LineConvex::Segment) {
            return false;
        }
    }
    //
    localLimits[0] = auxi_convexGrid::convertExtremitySegment(dbLocalLimits[0], drift, dx[DIM - 1]);
    localLimits[1] = auxi_convexGrid::convertExtremitySegment(dbLocalLimits[1], drift, dx[DIM - 1]);
    //
    if (voxelType == VoxelType::Composite) {
        localLimits[0]--;
        localLimits[1]++;
    }
    localLimits[0] = max(localLimits[0], gridLimits[DIM - 1][0]);
    localLimits[1] = min(localLimits[1], gridLimits[DIM - 1][1]);
    return (localLimits[0] < localLimits[1]);
}

inline void vox::auxi::auxi_convexGrid::verifySliceInstruction(vector<vox::auxi::SliceInstruction<long>> sliceInstruction) {
    bool problem = false;
    if (sliceInstruction.size() == 0) {
        return;
    }
    else {
        stable_sort(sliceInstruction.begin(), sliceInstruction.end(), [](auto& si1, auto& si2) {
            return si1.limits[0] < si2.limits[0];
            });
        long lmax = sliceInstruction[0].limits[1];
        for (size_t i = 1; i < sliceInstruction.size();i++) {
            if (not(sliceInstruction[i].limits[0] == lmax) or not(sliceInstruction[i].limits[0] < sliceInstruction[i].limits[1])) {
                problem = true;
            }
            lmax = sliceInstruction[i].limits[1];
        }
    }
    if (problem) {
        for (const auto& si : sliceInstruction) {
            cerr << si.limits[0] << " " << si.limits[1] << endl;
        }
        throw runtime_error("Big Problem");
    }
}

template<unsigned short DIM>
template<class C>
inline void auxi::ConvexGrid<DIM>::layerSlice(const C& singlePolyhedron, bool& firstPhaseNOTFound,
    const Point<DIM>& center, double distanceToSurface,
    array<long, 2>& limitFormerInclusion, size_t index, VoxelType voxelType, PhaseType phase) {
    const auto& ij = getAllFirstIndices()[index];
    array<long, 2> limitCurrentInclusion{};
    if (auxi_convexGrid::getLimits<DIM>(singlePolyhedron, center, ij, dx, limitCurrentInclusion, gridLimits, distanceToSurface * halfDiagonal, voxelType)) {
        array<long, 2> limitOfSlice;
        if (firstPhaseNOTFound) { // First inclusion. Here, the whole segment has to be a slice.
            limitOfSlice = limitCurrentInclusion;
            sliceInstructions[index] = { SliceInstruction<long>(limitOfSlice, voxelType, phase) };
            firstPhaseNOTFound = false;
            // update limitFomerInclusion
            limitFormerInclusion = limitCurrentInclusion;
        }
        else { // Other inclusions. Here, we have to add two slices. | slice1 | ... | slice2 |
            limitOfSlice[0] = limitCurrentInclusion[0];
            limitOfSlice[1] = limitFormerInclusion[0];
            if (limitOfSlice[0] < limitOfSlice[1]) sliceInstructions[index].push_back(SliceInstruction<long>(limitOfSlice, voxelType, phase));
            limitOfSlice[0] = limitFormerInclusion[1];
            limitOfSlice[1] = limitCurrentInclusion[1];
            if (limitOfSlice[0] < limitOfSlice[1]) sliceInstructions[index].push_back(SliceInstruction<long>(limitOfSlice, voxelType, phase));
            // update limitFomerInclusion
            limitFormerInclusion[0] = min(limitCurrentInclusion[0], limitFormerInclusion[0]);
            limitFormerInclusion[1] = max(limitCurrentInclusion[1], limitFormerInclusion[1]);
        }
    }
}

template<unsigned short DIM>
template<VoxelRule voxelRule, class C>
inline void auxi::ConvexGrid<DIM>::compute(const C& smallShape) {
    static_assert(voxelRule == VoxelRule::Center or voxelRule == VoxelRule::Average);
    //
    const auto& innerFaces = smallShape.getInnerInclusions();
    size_t indexMax = getAllFirstIndices().size();
    for (size_t index = 0; index < indexMax; index++) {
        sliceInstructions[index] = {};
        array<long, 2> limitFormerInclusion{};
        ///
        bool firstPhaseNOTFound = true;
        for (long iPhase = innerFaces.size() - 1; iPhase >= 0; iPhase--) {
            const auto& singlePolyhedron = innerFaces[iPhase];
            /// first case : Center
            if constexpr (voxelRule == VoxelRule::Center) {
                layerSlice(singlePolyhedron, firstPhaseNOTFound, smallShape.center, 0., limitFormerInclusion, index,
                    VoxelType::MonoPhase, smallShape.getPhaseForVoxellation(iPhase));
                // warning : not nice
                if constexpr (std::is_same<C, smallShape::SpheroPolyhedronInc<DIM>>::value) {
                    layerSlice(singlePolyhedron, firstPhaseNOTFound, smallShape.center, 1., limitFormerInclusion, index,
                        VoxelType::Composite, smallShape.getPhaseForVoxellation(iPhase));
                }
                // second case : Average
            }
            else if constexpr (voxelRule == VoxelRule::Average) {
                layerSlice(singlePolyhedron, firstPhaseNOTFound, smallShape.center, -1., limitFormerInclusion, index,
                    VoxelType::MonoPhase, smallShape.getPhaseForVoxellation(iPhase));
                layerSlice(singlePolyhedron, firstPhaseNOTFound, smallShape.center, 1., limitFormerInclusion, index,
                    VoxelType::Composite, smallShape.getPhaseForVoxellation(iPhase));
            }
        }
        // vox::auxi::auxi_convexGrid::verifySliceInstruction(sliceInstructions[index]); | FOR DEBUG debug
    }
}

} // namespace vox
} // namespace merope


#endif /* GRID_CONVEXGRID_IXX_ */
