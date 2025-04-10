//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"

#include "../../../AlgoPacking/include/SphereContainer.hxx"
#include "../../../AlgoPacking/include/ArrayDimensions.hxx"

#include "../SingleVoxel/SingleVoxel_Headers.hxx"
#include "../MesoStructure/InsideTorus.hxx"
#include "../MicroInclusion/MicroInclusion.hxx"


namespace merope {
namespace vox {

//! transforms nbNodes, L into dx
template<unsigned short DIM>
array<double, DIM> get_dx_from(array<size_t, DIM> nbNodes, array<double, DIM> L);
//! transforms dx, nbNodes into L
template<unsigned short DIM>
array<double, DIM> get_L_from(array<size_t, DIM> nbNodes, array<double, DIM> dx);

template<unsigned short DIM>
//! intermediate class storing dx
class With_dx {
public:
    //! constructor
    With_dx(array<double, DIM> dx_) : dx{ dx_ }, inverse_dx{}, halfDiagVoxel{ 0. } {
        setDx(dx);
    }
    //! getter
    const array<double, DIM>& getDx() const { return dx; }
    //! getter
    const array<double, DIM>& getInverseDx() const { return inverse_dx; }
    //! getter
    double getHalfDiagVoxel() const { return halfDiagVoxel; }

private:
    //! discretization step
    array<double, DIM> dx;
    //! inverse of discretization step
    array<double, DIM> inverse_dx;
    //! half diagonal of a voxel
    double halfDiagVoxel;
    //! \return the length of the half diagonal of a voxel
    double computeLengthHalfDiagVoxel() const { return 0.5 * sqrt(geomTools::normeCarre<DIM>(this->getDx())); }
protected:
    //! setter
    void setDx(const array<double, DIM>& dx_);
};

template<unsigned short DIM>
class GridParameters : public InsideTorus<DIM>, public SubArrayDimensions<DIM>, public With_dx<DIM> {
    //! the subgrid extracts a grid nMin<= ijk < nMax from the grid 0 <= ijk < nbNodes
public:
    GridParameters(array<size_t, DIM> nbNodes_, array<double, DIM> dx_);
    //! \warning : reset nMin and nMax
    void set_nbNodes_L(const array<size_t, DIM>& nbNodes_, const array<double, DIM>& L_);
private:
    using InsideTorus<DIM>::setLength;  // avoid setting L without setting dx
};

//! @return : a subgrid=grid of given dimensions
//! @param nbNodes : total number of voxels
//! @param L : lengths of the domain
template<unsigned short DIM>
GridParameters<DIM> create_grid_parameters_N_L(array<size_t, DIM> nbNodes, array<double, DIM> L);

//! @return : a subgrid of given dimensions
//! @param nbNodes : total number of voxels
//! @param nMin : minimal index of the subgrid.
//! @param nMax : maxmal index of the subgrid.
//! @param L : lengths of the domain
template<unsigned short DIM>
GridParameters<DIM> create_grid_parameters_N_L(array<size_t, DIM> nbNodes, array<size_t, DIM> nMin, array<size_t, DIM> nMax, array<double, DIM> L);


namespace auxi {
//! return the limits of the discrete subgrid containing the cuboid inside grid
//! \param cuboid : real space cuboid
//! \param preGrid : space discretisation
template<unsigned short DIM>
array<array<long, 2>, DIM> computeGridLimits(const Cuboid<DIM>& cuboid, const auto& preGrid);
//! return the intersection of part of the gridLimits that, once projected, is contained in the subgrid
//! \param torusGridLimits : gridLimits in the full grid
//! \param preSubGrid : subgrid parametrization
template<unsigned short DIM>
vector<array<array<long, 2>, DIM>> intersectGridLimits(const array<array<long, 2>, DIM>& gridLimits_, const GridParameters<DIM>& preSubGrid);
//! recovers a sequence of segments such that their projections by periodicity lie in nbMin <= i < nMax
vector<array<long, 2>> auxi_intersectGridLimits(const array<long, 2>& gridLimits, long nbNodes, long nMin, long nMax);
}  // namespace auxi

}  // namespace vox
}  // namespace merope

#include "../Grid/AreCompatible.hxx"
#include "../Grid/PreGrid.ixx"


