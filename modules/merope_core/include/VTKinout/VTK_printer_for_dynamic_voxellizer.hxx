//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "../VTKinout/VTK_adapter.hxx"
#include "../Voxellation/DynamicVoxellizer.hxx"

namespace merope {
namespace vtk_adapter {

template<unsigned short DIM>
class vtk_printer {
public:
    inline void printVTK(const vox::voxellizer::GridRepresentation<DIM>& gridRepresentation,
        string fileVTK,
        string nameValue = "MaterialId") {
        auto pointer_to_cartesian_grid = gridRepresentation.template get_if<PhaseType>();
        if (pointer_to_cartesian_grid) {
            vtk_adapter::printVTK<unsigned short>(*pointer_to_cartesian_grid, fileVTK, nameValue);
            return;
        }
        auto pointer_to_cartesian_grid_2 = gridRepresentation.template get_if<double>();
        if (pointer_to_cartesian_grid_2) {
            vtk_adapter::printVTK<double>(*pointer_to_cartesian_grid_2, fileVTK, nameValue);
            return;
        }
        gridRepresentation.raise_error_internal_type_message(__PRETTY_FUNCTION__);
    }

    inline void printVTK_segmented(const vox::voxellizer::GridRepresentation<DIM>& gridRepresentation,
        string fileVTK, string fileCoeff,
        string nameValue = "MaterialId") {
        auto pointer_to_cartesian_grid = gridRepresentation.template get_if<double>();
        if (pointer_to_cartesian_grid) {
            vtk_adapter::printVTK_segmented<unsigned short>(*pointer_to_cartesian_grid,
                fileVTK, fileCoeff, nameValue);
        } else {
            gridRepresentation.raise_error_internal_type_message(__PRETTY_FUNCTION__);
        }
    }

    inline void printVTK_removeUnusedPhase(const vox::voxellizer::GridRepresentation<DIM>& gridRepresentation,
        string fileVTK, string fileCoeff,
        string nameValue = "MaterialId") {
        auto pointer_to_cartesian_grid = gridRepresentation.template get_if<PhaseType>();
        if (pointer_to_cartesian_grid) {
            vtk_adapter::printVTK_removeUnusedPhase<unsigned short>(*pointer_to_cartesian_grid,
                fileVTK, fileCoeff, nameValue);
        } else {
            gridRepresentation.raise_error_internal_type_message(__PRETTY_FUNCTION__);
        }
    }
};


}  // namespace  vtk_adapter
}  // namespace  merope