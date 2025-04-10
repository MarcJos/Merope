//! Copyright : see license.txt
//!
//! \briefSpecial class for providing clever access to corners, knowing the relative position of the sphere
//! w.r.t. the cube, in order to see whether it intersects it or not.
//!
//! \fixme : try to make it less messy.

#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "../../GenericTools/Loops.hxx"

#include "../../Geometry/include/AmbiantSpace.hxx"

namespace sac_de_billes {

class Path {
    Path(Path&&) = delete;
    Path(const Path&) = delete;
    Path& operator=(Path&&) = delete;
    Path& operator=(const Path&) = delete;
    Path();

public:
    static const Path& get();
    // get all the path for neighboring cubes
    template<unsigned DIM, size_t NB_NGHB>
    const array<array<int, DIM>, nbSubcubes<DIM>(NB_NGHB)>& pathForNeighbors() const;
    //! path for corners
    const vector<vector<array<unsigned short, 3>>> pathForCorners;
    //! efficient storage of corner information by a single vector entry
    const vector<vector<array<unsigned short, 3>>> pathFor27Corners;
    //! 27 possibilities of relative position wrt a cube
    const array<array<int, 3>, 27> CORNER_27_DIRS;
    //! efficient access to corners stored in pathFor27Corners
    const vector<array<unsigned short, 3>>& myCornersEfficient(size_t indexPosRelative) const;
    //! backward transformation from spatial corner localization to index referring to pathFor27Corners.
    template<unsigned short DIM>
    size_t fromCorner2Index(array<int, DIM> posRelative) const;

private:
    static constexpr size_t NB_NEIB_MAX = 5;
    static constexpr array<int, 7> IX = { 0, 1, -1, 2, -2, 3, -3 };
    //! precomputed pathForNeighbors
    const array<array<int, 3>, nbSubcubes<3>(5)> pathNghb_3D_5;
    const array<array<int, 3>, nbSubcubes<3>(7)> pathNghb_3D_7;
    const array<array<int, 2>, nbSubcubes<2>(5)> pathNghb_2D_5;
    const array<array<int, 2>, nbSubcubes<2>(7)> pathNghb_2D_7;

    vector<array<unsigned short, 3>> myCorners(array<int, 3> posRelative) const;

    vector<vector<array<unsigned short, 3>>> buildPathForCoin();
    //! efficient way to test spheres, from the nearest voxel to the farthest
    template<unsigned short DIM, size_t NB_NGHB>
    array<array<int, DIM>, nbSubcubes<DIM>(NB_NGHB)> createParcours();
    //! 27 possibilities of relative position wrt a cube
    //! for each index : -1 < 0 = same plane < 1
    array<array<int, 3>, 27> create27Directions() const;
    //! Creates the vector containing the relevant corners
    //! that shall be tested by a sphere positionned with a relative position relPose
    vector<vector<array<unsigned short, 3>>> build27Corners();
};

namespace path {
template<unsigned short DIM>
class TabCorner {
    TabCorner(TabCorner&&) = delete;
    TabCorner(const TabCorner&) = delete;
    TabCorner& operator=(TabCorner&&) = delete;
    TabCorner& operator=(const TabCorner&) = delete;
    TabCorner();
public:
    static const TabCorner<DIM>& get();
    // opposite corners in order to reach get more easily out of the sphere
    const vector<array<unsigned short, DIM>>& getTab() const { return storage_tab_corner; }
private:
    vector<array<unsigned short, DIM>> storage_tab_corner;
};

}  // namespace  path

}  // namespace sac_de_billes

#include "Path.ixx"


