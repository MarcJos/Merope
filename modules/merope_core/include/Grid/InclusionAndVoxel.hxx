//! Copyright : see license.txt
//!
//! \brief
//!
//
#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Grid/GridTypes.hxx"
#include "../Grid/VoxelRule.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

namespace vox {
namespace inclusionAndVoxel {

//! simple method for filling a voxel with given data
//! \param voxelData : data of the voxel to be updated
//! \param phases2Include : data to feed in the voxelData (either replace or push_back, depending on the composite::OutputFormat)
template<class VOXEL_TYPE>
void fillVoxel(VOXEL_TYPE& voxelData, const VOXEL_TYPE& phases2Include);
//! simple method for filling a voxel with given data from an inclusion. Return if the voxel has been filled
//! \param microInclusion : inclusion considered
//! \param voxelData : data of the voxel to be updated
//! \param dx : dimensions of a voxel
//! \param halfDiagVoxel : length of the half-diagonal of the voxel
template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
bool fillVoxel(const INCLUSION& microInclusions, const DiscPoint<DIM>& indexVoxel, const Point<DIM>& dx, const double& halfDiagVoxel, VOXEL_TYPE& voxelData);
//! \return whether the microInclusion intersects the voxel. If yes, the data describing the intersection is inside phases2Include
//! \param phases2Include : data describing the intersection of the voxel and the microInclusion (if applicable)
//! \param indexVoxel : 3-D position of the voxel
//! \param microInclusion : microInclusion
//! \param dx : dimensions of the voxel
//! \param halfDiagVoxel : stores the lenght of the half-diagonal of the voxel
template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
inline bool phasesInsideVoxel(const INCLUSION& microInclusion, const DiscPoint<DIM>& indexVoxel, const Point<DIM>& dx, const double& halfDiagVoxel,
        VOXEL_TYPE& phases2Include);
//! \return the data consisting of the intersection of a shape with layers, and a voxel. If no intersection, return false.
//! \param smallShape : a big inclusion
//! \param centerPoly_to_origVoxel : a vector from the center of the smallShape to the origin of the voxel
//! \param dx : dimensions of the voxel
//! \param halfDiagVoxel : (measure of the) half-diagonal of the voxel
template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
inline VOXEL_TYPE computeAllFracVol(const INCLUSION& inclusion, const Point<DIM>& centerPoly_to_origVoxel, const Point<DIM>& dx, const double& halfDiagVoxel);
// tolerance above which the result is considered close to 1
constexpr double TABCOEFFS_TOL = 0.999;

}  // namespace inclusionAndVoxel
}  // namespace vox
}  // namespace merope


#include "../Grid/InclusionAndVoxel.ixx"


