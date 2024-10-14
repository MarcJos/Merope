//! Copyright : see license.txt
//!
//! \brief
//!
//!

#pragma once


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../Grid/Grid.hxx"
#include "../Grid/GridTypes.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"
#include "../Physics/Homogenization.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {

// for linearizing a integer multi-index (i[0], i[1], i[2]) into an grid of dimensions (dims[0], dims[1], dims[2]) into a single index
enum class IndexConvention {
    AMITEX,  // convention used by AMITEX:  i[0] + (i[1] + i[2] * dims[1]) * dims[0]
    Merope    // convention used by Merope or Mérope  i[2] + (i[1] + i[0] * dims[1]) * dims[2]
};

class Grid_VER : public Grid {
public:
    //! Default constructor
    Grid_VER();
    //! Grid constructor by reading a periodical continuous geometry
    //! \param geometry the periodical discrete geometry
    Grid_VER(VTKRead& geometry);
    //! Grid constructor by reading a TIFF File
    //! \param tf: TIFF File
    //! \param dx: Grid step
    Grid_VER(TIFFRead& tf, double);

    array<size_t, 3> nbNodes;
    vector<double> L;
    array<double, 3> dx;
    array<double, 3> inverse_dx;

    //! fixme
    template<unsigned short DIM>
    void fromGridPhase(const vox::GridPhase<DIM>&, vector<double>& coefficients);
    //! fixme
    void removeUnusedPhase(vector<double>& coefficients);

    //! symmetrise the grid in all the directions
    template<unsigned short DIM>
    void symmetrize(array<size_t, DIM> symDirections);
    //! fixme
    template<unsigned short DIM>
    void set_Nb(array<size_t, DIM> nb);
    //! set dimension and lengths of the torus
    //! \param L_ : lengths of the torus
    template<unsigned short DIM>
    void set_L(array<double, DIM> L_);

    //! fixme
    template<unsigned short DIM, IndexConvention indexConvention = IndexConvention::Merope>
    size_t get_linear_index(const array<size_t, DIM>& ijk) const;
    //! fixme
    size_t get_linear_index(size_t i, size_t j, size_t k) const;
    //! fixme
    template<unsigned short DIM>
    array<size_t, DIM> get_coord_index(size_t i) const;
    //! @tparam indexConvention : way of flattening 3D indices into a single index
    //! @return : a vector of {zones}, each zone containing all the indices of voxels of the same type
    template<unsigned short DIM, IndexConvention indexConvention>
    vector<vector<size_t>> get_phase_list() const;

private:
    //! fixme
    void set_Voxel(const array<size_t, 3>& ijk, size_t indexPhase);
    //! fixme
    void set_Voxel(size_t i, size_t j, size_t k, size_t indexPhase);
    //! fixme
    void set_Voxel(size_t linear_index, size_t indexPhase);
    //! fixme
    void set_NbPhase(size_t s);

    //! Common initialisation
    void initCommon();
};

namespace vox {
//! symmetrize a grid
template<unsigned short DIM>
void symmetrize(Grid_VER& myGrid, array<size_t, DIM> symDirections);
//! symmetrize a grid
template<unsigned short DIM>
void symmetrize(string inputFileName, string outputFileName, array<size_t, DIM> symDirections);

namespace aux {
//! see Grid_VER::symmetrize
size_t symmetriseAuxi(size_t N1, size_t i, size_t N2);
}  // namespace aux

}  // namespace vox
}  // namespace merope


#include "../Grid/Grid_VER.ixx"


