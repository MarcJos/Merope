//! Copyright : see license.txt
//!
//! \brief

#ifndef VTKINOUT_VTK_ADAPTER_HXX
#define VTKINOUT_VTK_ADAPTER_HXX


#include "../MeropeNamespace.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "VTKStream.hxx"

namespace merope {
template<unsigned short DIM, class VOXEL_TYPE>
void printVTK(string fileVTK, const vox::CartesianGrid<DIM, VOXEL_TYPE>& cartesianGrid,
    string nameValue, string typeValue);
} // namespace  merope

#include "VTK_adapter.ixx"

#endif // VTKINOUT_VTK_ADAPTER_HXX