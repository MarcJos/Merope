//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef PREVOXELLATION_HXX_
#define PREVOXELLATION_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Geometry/GeomTools.hxx"
#include "../Grid/Grid_VER.hxx"
#include "../Grid/PreGrid.hxx"
#include "../Grid/VoxelRule.hxx"
#include "../MultiInclusions/LaguerreTess.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"
#include "../MultiInclusions/SphereInclusions.hxx"
#include "../MesoStructure/Structure.hxx"


#include "../MeropeNamespace.hxx"


namespace merope {
namespace vox {

template<unsigned short DIM>
class PreVoxellation: public PreSubGrid<DIM>, public WithVoxelRule {
public:
    PreVoxellation(array<double, DIM> L, array<size_t, DIM> nbNodes_ = create_array<DIM, size_t>(1), VoxelRule voxelRule_ =
        VoxelRule::Center):
        PreSubGrid<DIM>(nbNodes_, vox::get_dx_from<DIM>(nbNodes_, L)), WithVoxelRule(voxelRule_), grid{} {
        this->setGridNL(nbNodes_, L);
        this->setSubGridIndices(create_array<DIM, size_t>(0), nbNodes_);
    }
    //! gets the internal grid
    const Grid_VER& getGrid() const { return grid; }

protected:
    //! stores the grid
    Grid_VER grid;
    //! sets grid characteristics
    void setGridNL(array<size_t, DIM> nbNodes, array<double, DIM> L);
};

} // namespace vox
} // namespace merope


#include "../Voxellation/PreVoxellation.ixx"

#endif /* PREVOXELLATION_HXX_ */
