//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../../AlgoPacking/src/StdHeaders.hxx"
#include "../MeropeNamespace.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"

#include "../Geometry/GeomTools.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../Grid/GridTypes.hxx"
#include "../Grid/PreGrid.hxx"
#include "../Grid/VoxelRule.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"

namespace merope {

namespace vox {
namespace auxi {

enum class VoxelType {
        MonoPhase, Composite
};

template<class INT_TYPE>
struct SliceInstruction {
        //! Instruction for filling a grid in a slice directed by the last coordinate.
        //! The line is directed by fixing the first coordinates (not recorded here)
        //! 2 possibilities:
        //! -> we know the phase, recorded into phase
        //! -> we do NOT know the phase, and then have to consider composite voxels
        SliceInstruction(array<INT_TYPE, 2> limits_, VoxelType voxelType_, PhaseType phaseId_) :
                limits(limits_), voxelType(voxelType_), phase(phaseId_) {}
        //! limits of the slice
        array<INT_TYPE, 2> limits;
        //! type of voxel : MonoPhase or Composite
        VoxelType voxelType;
        //! phase
        PhaseType phase;
        //! print
        void print(std::ostream& ost) {
                ost << limits[0] << " : " << limits[1] << " ; "
                        << static_cast<int>(voxelType) << " ; " << phase;
        }
};


template<unsigned short DIM>
class ConvexGrid {
        //! class for implementing efficient voxellation of convex shapes.
        //! The principle is the following, when voxellizing (in 3D), we inspect the voxels{i,j,k} and decide is they are inside a given shape
        //! Since the shape is convex, this can be turnt into
        //! Finding {k_min(i,j), k_max(i,j)} such that any voxel {i,j,k} for k \in {k_min(i,j), k_max(i,j)} is inside the shape
        //! The whole purpose of this class is to find and expose the indices {k_min, k_max}
public:
        ConvexGrid(const array<array<long, 2>, DIM>& gridLimits, const array<double, DIM>& dx);
        //! \return compute the vector limits == list of {k_min,kmax}
        //! \param C : a SmallShape (intrusive on which shape it is)
        template<VoxelRule voxelRule, class C>
        void compute(const C& smallShape);
        //! getter
        size_t getNbFirstIndices() const { return allFirstIndices.size(); }
        //! \return a DiscPoint<DIM> {i,j,k}, with i,j corresponding to allFirstIndices[index]
        //! \param index : index of the {i,j} firstindices in the ConvexGrid
        const DiscPoint<DIM>& getIndexAllCoordinates(size_t index) const { return allFirstIndices[index]; }
        //! \return the vector of whole slice instructions corresponding to a fixed index
        //! \param index : index of the {i,j} firstindices in the ConvexGrid
        vector<SliceInstruction<long>> getSliceInstructions(size_t index) const { return sliceInstructions[index]; }
private:
        //! helper function. It adds new slices to the sliceInstructions, corresponding to the surface
        //! \param singleSurface : encodes the surface of the inclusion
        //! \param firstPhaseNOTFound : tells whether the first inclusion intersecting the voxel line
        //! \param distanceToSurface : modifies the distance to the surface
        //! \param limitFormerSurface : stores the slice indices of the former surface
        //! \param index : index of the voxel line
        //! \param voxelType : type of the voxel (composite or monophase)
        //! \param phase : identifies the phase
        template<class C>
        void layerSlice(const C& singleSurface, bool& firstPhaseNOTFound,
                const Point<DIM>& center, double distanceToSurface,
                array<long, 2>& limitFormerSurface, size_t index, VoxelType voxelType, PhaseType phase);
        //! getter
        const vector<DiscPoint<DIM>>& getAllFirstIndices() const { return allFirstIndices; }

        //! space discretization
        array<double, DIM> dx;
        //! norm of the halfDiagonal
        double halfDiagonal;
        //! contains the limits of the grid, in the form for {{lmin[0], lmax[0]}, {lmin[1],lmax[1]}, ...}
        array<array<long, 2>, DIM> gridLimits;
        //! contains the number of voxels in each direction i in {0,DIM-1}
        array<size_t, DIM> gridSize;
        //! contains all the pairs {i,j} such that there exists a voxel {i,j,k} in the grid
        const vector<DiscPoint<DIM>> allFirstIndices;
        //! vector containing all the limits k_min< l < k_max corresponding to a certain phase. Each limits[index] correspond to a row {i,j} = allFirstIndices[index]
        vector<vector<SliceInstruction<long>>> sliceInstructions;
};

namespace auxi_convexGrid {
//! \return a long i with the property that the segment [i(convertExtremitySegment[x]), i(convertExtremitySegment[y])) * dx in included into [x+drift, y+drift)
//! \param realSegment : the initial segment [a_0, a_1]
//! \param drift : the drift, ie we should consider [realSegment[0]+drift, realSegment[1]+drift]
//! \param dx : the dx step in the desired direction
inline long convertExtremitySegment(double extremitySegment, double drift, double dx) { return static_cast<long>(floor((extremitySegment + drift) / dx)) + 1; }


// WARNING : this section, although correct, is not clear at all, and may provoke bugs. Everything should be clarified and commented in a factorized way.

//! \return true if the smallshape intersects the voxel line x = x1x2[0] and y = x1x2[1] (in 3D)
//! \param smallShape : vector<Halspace> describing a polyhedron or a sphere
//! \param x1x2 : 1st and 2nd coordinates of a line (in 3D)
//! \param ij : 1st and 2nd coordinates of a line of voxels (in 3D)
//! \param dx : dimensions of each voxel in the line
//! \param currentLimits : on line (x,y)=x1x2, 
//!                             (i) [k_min, k_max) *dx[DIM-1] inside the smallShape
//!                             (ii) (k_min-1) * dx[DIM-1] and kmax * dx[DIM-1] outside the shape
//!                             if distance > 0 (ii) is guaranteed
//!                             if distance < 0 (i) is guaranteed
//!                             WARNING : in general, (i) and (ii) are not simulatenously guaranteed
//! \param gridLimits: stores the vector {{i, i+1}, {j,j+1}, {k_min, k_max}} that describes the intersection
//! \param distance : distance to make the shape grow othogonally to its surface. (Negative = inside the shape, positive = outside the shape).
template<unsigned short DIM, class C>
bool getLimits(const C& smallShape, const Point<DIM>& center, const DiscPoint<DIM>& ij,
        const array<double, DIM>& dx, array<long, 2>& currentLimits,
        const array<array<long, 2>, DIM>& gridLimits, double distance, VoxelType voxelType);


//! compute the intersection of a voxel line i = ij[0] and j = ij[1] (in 3D) with a convex inclusion
//! \param dbLocalLimits : store the limits. The new limits should be inside the given limits
//! \param inclusion: inclusion to interset
//! \param x1x2 : the line is directed in the 3rd canonical direction, and goes through {x1x2[0], x1x2[1], 0}
//! \param distance : distance to make the shape grow othogonally to its surface. 
//! = 0 : guarantees the intersection of the line with the shape is exact
//! < 0 : guarantees the computed line is inside the set of points inside the shape at a distance at least (param)distance of the shape
//! > 0 : guarantees the intersection of the shape with the computed line contains the set of points of the line either inside the shape or outside it as a distance at most (param)distance of the shape
template<unsigned short DIM, class Inclusion>
geomTools::Intersection_LineConvex getLimits_innerFunction(array<double, 2>& dbLocalLimits, const Inclusion& inclusion,
        const Point<DIM - 1>& x1x2, double distance);

//! fixme
void verifySliceInstruction(vector<vox::auxi::SliceInstruction<long>>);
}  // namespace auxi_convexGrid

}  // namespace auxi
}  // namespace vox
}  // namespace merope

#include "../Grid/ConvexGrid.ixx"


