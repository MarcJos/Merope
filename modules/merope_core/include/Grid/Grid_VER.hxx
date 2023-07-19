//! Copyright : see license.txt
//!
//! \brief 
//!
//!

#ifndef GRID_GRID_VER_HXX_
#define GRID_GRID_VER_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"
#include "../Grid/CartesianGrid.hxx"
#include "../Grid/Grid.hxx"
#include "../Grid/GridTypes.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"
#include "../Physics/Homogenization.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM>
class Voxellation;
}


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
    //! Print the internal variables
    //! \param os: Output stream
    Grid_VER(const vox::Voxellation<3>& voxellisation);

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
    template<unsigned short DIM>
    size_t get_linear_index(const array<size_t, DIM>& ijk) const;
    //! fixme
    size_t get_linear_index(size_t i, size_t j, size_t k) const;
    //! fixme
    template<unsigned short DIM>
    array<size_t, DIM> get_coord_index(size_t i) const;

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
} // namespace aux

} // namespace vox
} // namespace merope


#include "../Grid/Grid_VER.ixx"

#endif /* GRID_GRID_VER_HXX_ */
