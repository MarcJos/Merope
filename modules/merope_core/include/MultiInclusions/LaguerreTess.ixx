//! Copyright : see license.txt
//!
//! \brief

#pragma once


#include "../../../Geometry/include/GeomTools.hxx"

#include "../Voronoi/VoroInterface.hxx"


namespace merope {

// LaguerreTess

template<unsigned short DIM>
inline LaguerreTess<DIM>::LaguerreTess(array<double, DIM> L_, vector<Sphere<DIM> > seeds_) :
    PolyInclusions<DIM>(), WithAspratio<DIM>(), seeds{ seeds_ } {
    this->setLength(L_);
}

template<unsigned short DIM>
inline void LaguerreTess<DIM>::computeTessels() {
    bool trivial_aspratio = true;
    for (int i = 0; i < DIM; i++) {
        trivial_aspratio = trivial_aspratio
            and (abs(1 - this->aspratio[i]) < geomTools::EPSILON);
    }
    if (trivial_aspratio)
        no_aspratio_computeTessels();
    else with_aspratio_computeTessels();
}

template<unsigned short DIM>
inline void LaguerreTess<DIM>::no_aspratio_computeTessels() {
    voroInterface::VoroInterface<DIM> voroInterface(this->tore.L, seeds);
    this->setInclusions(voroInterface.getMicroInclusions());
    for (size_t i = 0; i < this->getMicroInclusions().size(); i++) {  // recovers the correct phase
        this->getMicroInclusions()[i].getPhaseGraphical(0) = seeds[i].phase;
    }
}

template<unsigned short DIM>
inline void LaguerreTess<DIM>::with_aspratio_computeTessels() {
    array<double, DIM> L_dilated = this->tore.L;
    vector < Sphere <DIM>> seeds_dilated = seeds;
    // transforms in the dilated coordinates
    linearTransform::point <DIM>(L_dilated, this->inverse_aspratio);
    for (auto& sph : seeds_dilated) {
        linearTransform::point <DIM>(sph.center, this->inverse_aspratio);
    }
    // transforms back in the original coordinates
    LaguerreTess <DIM> res(L_dilated, seeds_dilated);
    res.computeTessels();
    this->setInclusions(res.getMicroInclusions());
    for (auto& tess : this->getMicroInclusions()) {
        tess.linearTransform(this->aspratio);
    }
}

}  // namespace merope



