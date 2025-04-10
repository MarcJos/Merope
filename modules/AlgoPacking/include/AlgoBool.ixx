//! Copyright : see license.txt
//!
//! \briefImplements random boolean sphere placement
//
#pragma once

#include "RandomShooter.hxx"

namespace sac_de_billes {

template<unsigned short DIM>
inline AlgoBool_Template<DIM>::AlgoBool_Template(
    AmbiantSpace::BigShape<DIM>* T_,
    algoRSA_aux::RadiusGenerator<DIM>* radiusGen_, unsigned seed_) :
    placedSpheres{ new vector<Sphere<DIM>> { } }, bigShape{ T_ }, radiusGen{
            radiusGen_ } {
    // Intialize the random generator, cf https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution
    mt19937 gen(seed_);  //Standard mersenne_twister_engine seeded with rd()
    randGenerator = gen;
}

template<unsigned short DIM>
inline map<string, string> AlgoBool_Template<DIM>::proceed(int) {
    bool oneMoreTime = true;
    while (oneMoreTime) {
        auto newPoint = pickPoint();
        auto newSphere = Sphere < DIM
        >(newPoint, radiusGen->getRadius(), radiusGen->getPhase());
        if (bigShape->isInside(newSphere)) {
            placedSpheres->push_back(newSphere);
            oneMoreTime = radiusGen->nextRadius();
        }
    }
    return map<string, string>{ { "Packed", "False" }};
}

template<unsigned short DIM>
inline Point<DIM> AlgoBool_Template<DIM>::pickPoint() {
    return randomShooter::pickInCuboid<DIM>(randGenerator, bigShape->L);
}

}  // namespace sac_de_billes


