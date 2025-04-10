//! Copyright : see license.txt
//!
//! \brief

#include "Point.hxx"

#pragma once

namespace merope {
namespace vox {
namespace auxi {
//! \return the origin of a voxel, in continuous coordinates designated by its discrete origin
//! \param discPt : discrete coordinates of the origini of the voxel
//! \param dx : lengths of the voxels
template<unsigned short DIM, class C>
Point<DIM> origin(const C& discPt, const Point<DIM>& dx);
//! \return the center of a voxel, in continuous coordinates
//! \param discPt : discrete coordinates of the origini of the voxel
//! \param dx : lengths of the voxels
template<unsigned short DIM>
Point<DIM> center(const DiscPoint<DIM>& discPt, const Point<DIM>& dx);
}  // namespace  auxi
}  // namespace  vox
}  // namespace  merope


#include "Voxel.ixx"