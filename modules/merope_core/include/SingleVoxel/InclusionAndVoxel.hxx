//! Copyright : see license.txt
//!
//! \brief
//!
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"
#include "../SingleVoxel/VoxelRule.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"


namespace merope {

namespace vox {
namespace inclusionAndVoxel {

//! simple method for filling a voxel with given data from an inclusion. Return if the voxel has been filled
//! \param microInclusion : inclusion considered
//! \param voxelData : data of the voxel to be updated
//! \param dx : dimensions of a voxel
//! \param halfDiagVoxel : length of the half-diagonal of the voxel
template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE, class VOXEL_POLICY>
bool fillVoxel(const INCLUSION& microInclusions, const DiscPoint<DIM>& indexVoxel, const Point<DIM>& dx, const double& halfDiagVoxel, VOXEL_TYPE& voxelData,
        const VOXEL_POLICY& voxelPolicy);
//! \return whether the microInclusion intersects the voxel. If yes, the data describing the intersection is inside phases2Include
//! \param phases2Include : data describing the intersection of the voxel and the microInclusion (if applicable)
//! \param indexVoxel : 3-D position of the voxel
//! \param microInclusion : microInclusion
//! \param dx : dimensions of the voxel
//! \param halfDiagVoxel : stores the lenght of the half-diagonal of the voxel
template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
bool phasesInsideVoxel(const INCLUSION& microInclusion, const DiscPoint<DIM>& indexVoxel, const Point<DIM>& dx, const double& halfDiagVoxel,
        VOXEL_TYPE& phases2Include);
//! \return the data consisting of the intersection of a shape with layers, and a voxel. If no intersection, return false.
//! \param inclusion : a big inclusion
//! \param centerPoly_to_origVoxel : a vector from the center of the smallShape to the origin of the voxel
//! \param dx : dimensions of the voxel
//! \param halfDiagVoxel : (measure of the) half-diagonal of the voxel
template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
VOXEL_TYPE computeAllFracVol(const INCLUSION& inclusion, const Point<DIM>& centerPoly_to_origVoxel, const Point<DIM>& dx, const double& halfDiagVoxel);

//! \return the data consisting of the intersection of a shape with layers, and a voxel. If no intersection, return false.
//! \param inclusion : a big inclusion
//! \param centerPoly_to_origVoxel : a vector from the center of the smallShape to the origin of the voxel
//! \param dx : dimensions of the voxel
template<unsigned short DIM, class INCLUSION, class VOXEL_TYPE>
VOXEL_TYPE computeSimpleGeometry(const INCLUSION& inclusion, const Point<DIM>& centerPoly_to_origVoxel, const Point<DIM>& dx);

//! @return : the list 
//! \param inclusion : a big inclusion
//! \param centerPoly_to_origVoxel : a vector from the center of the smallShape to the origin of the voxel
//! \param dx : dimensions of the voxel
//! \param layer_i : index of the layer
template<unsigned short DIM, class INCLUSION>
tuple<bool, vector<HalfSpace<DIM>>> get_list_tangentPlanes(const INCLUSION& inclusion, const Point<DIM>& centerPoly_to_origVoxel, const Point<DIM>& dx, size_t layer_i);


// tolerance above which the result is considered close to 1
constexpr double TABCOEFFS_TOL = 0.999;

}  // namespace inclusionAndVoxel
}  // namespace vox
}  // namespace merope


#include "InclusionAndVoxel.ixx"


