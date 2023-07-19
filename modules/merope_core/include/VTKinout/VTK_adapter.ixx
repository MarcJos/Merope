//! Copyright : see license.txt
//!
//! \brief

#ifndef VTKINOUT_VTK_ADAPTER_IXX
#define VTKINOUT_VTK_ADAPTER_IXX


namespace merope {
template<unsigned short DIM, class VOXEL_TYPE>
void printVTK(string fileVTK, const vox::CartesianGrid<DIM, VOXEL_TYPE>& cartesianGrid,
    string nameValue, string typeValue) {
    VTKstream vtkstream(fileVTK.c_str());
    // header
    auto n = cartesianGrid.getNbNodeSubGrid();
    auto dx = cartesianGrid.getDx();
    if constexpr (DIM == 2) {
        vtkstream.STRUCTURED_POINTS(n[0], n[1], dx[0], dx[1]);
    } else if constexpr (DIM == 3) {
        vtkstream.STRUCTURED_POINTS(n[0], n[1], n[2], dx[0], dx[1], dx[2]);
    }
    size_t ng = auxi_function::productOf<size_t>(n);
    vtkstream.setCELL(ng);
    //
    vtkstream << "SCALARS " << nameValue << " " << typeValue << endl;
    vtkstream << "LOOKUP_TABLE default" << endl;
    vtkstream.writeVector(cartesianGrid, n);
}

} // namespace  merope

#endif // VTKINOUT_VTK_ADAPTER_IXX