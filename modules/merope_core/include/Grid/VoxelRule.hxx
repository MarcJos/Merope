//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../MeropeNamespace.hxx"

namespace merope {
namespace vox {

//! parameters for computing the voxellation
enum class VoxelRule {
    Laminate,   // compute the percentage of each phase inside the voxel + normal
    Average,    // compute the percentage of each phase inside the voxel
    Center      // evaluate the phase at the very center of the voxel
};

}  // namespace vox
}  // namespace merope


