//! Copyright : see license.txt
//!
//! \brief

#include "../../Geometry/include/VolumeInCube.hxx"
#include "../../Geometry/include/GeomTools.hxx"
#include "voro++.hh"

namespace merope {
namespace geomTools {

double volume_in_cube_3(const vector<HalfSpace<3>>& list_hf) {
    if (list_hf.size() == 0) {
        return 1.;
    } else if (list_hf.size() == 1) {
        return fracVolIntersection<3>(list_hf[0]);
    } else { // (list_hf.size() > 1)
        // Initialize the Voronoi cell to be a cube of side length 2, centered
        // on the origin
        voro::voronoicell poly;
        poly.init(0., 1., 0., 1, 0., 1.);
        //
        bool deleted_cell = false; // if cell is deleted, volume = 0 and stop
        for (const auto& hf : list_hf) {
            deleted_cell = not(poly.plane(hf.vec()[0], hf.vec()[1], hf.vec()[2], 2 * hf.c()));
            if (deleted_cell) {
                break;
            }
        }
        double volume = deleted_cell ? 0 : poly.volume();
        return volume;
    }
}

}  // namespace  geomTools
}  // namespace  merope
