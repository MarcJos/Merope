//! Copyright : see license.txt
//!
//! \brief

#pragma once

#include "GeomTools.hxx"
#include "GeomTools_1.hxx"


namespace merope {
namespace geomTools {

double volume_in_cube_3(const vector<HalfSpace<3>>& list_hf);

//! @return the volume of the intersection of a list of half-spaces and a given cube
//! @param list_hf : list of half-spaces
template<unsigned short DIM>
double volume_in_cube(const vector<HalfSpace<DIM>>& list_hf) {
    double volume = 0.;
    if constexpr (DIM == 2) {
        vector<HalfSpace<3>> list_hf_3D = {};
        for (size_t i = 0; i < list_hf.size(); i++) {
            auto vec = list_hf[i].vec();
            list_hf_3D.push_back(HalfSpace<3>({ vec[0], vec[1], 0 }, list_hf[i].c()));
        }
        return volume_in_cube_3(list_hf_3D);
    } else if constexpr (DIM == 3) {
        volume = volume_in_cube_3(list_hf);
    } else {
        Merope_static_error(HalfSpace<DIM>, "Incorrect dimension");
    }
    return volume;
}


}  // namespace  geomTools
}  // namespace  merope
