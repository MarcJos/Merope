//! Copyright : see license.txt
//!
//! \brief Implements random boolean sphere placement
//
#ifndef ALGOBOOL_IXX_
#define ALGOBOOL_IXX_

#include "RandomShooter.hxx"

namespace sac_de_billes {

template<unsigned short DIM>
inline AlgoBool_Template<DIM>::AlgoBool_Template(
    AmbiantSpace::BigShape<DIM>* T_,
    algoRSA_aux::RadiusGenerator<DIM>* radiusGen_, unsigned seed_):
    placedSpheres{ new vector<Sphere<DIM>> { } }, bigShape{ T_ }, radiusGen{
            radiusGen_ } {
    // Intialize the random generator, cf https://en.cppreference.com/w/cpp/numeric/random/uniform_int_distribution
    mt19937 gen(seed_); //Standard mersenne_twister_engine seeded with rd()
    randGenerator = gen;
    randomReal = uniform_real_distribution<>(0., 1.);
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
    return randomShooter::pickInCuboid<DIM>(randGenerator, randomReal, bigShape->L);
}

} // namespace sac_de_billes

#endif /* ALGOBOOL_IXX_ */
