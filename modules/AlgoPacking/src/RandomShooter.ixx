//! Copyright : see license.txt
//!
//! \brief
//
#ifndef RANDOMSHOOTER_IXX_
#define RANDOMSHOOTER_IXX_

namespace sac_de_billes {
using namespace std;

namespace randomShooter {

template<unsigned short DIM>
inline Point<DIM> pickInCuboid(mt19937& randGenerator,
    uniform_real_distribution<>& randomReal,
    const Point<DIM>& length) {
    Point <DIM> point{};
    for (size_t i = 0; i < DIM; i++) {
        point[i] = randomReal(randGenerator) * length[i];
    }
    return point;
}

template<unsigned short DIM>
inline Point<DIM> pickOnSphere(mt19937& randGenerator,
    uniform_real_distribution<>& randomReal) {
    Point<DIM> length = create_array<DIM>(1.);
    Point<DIM> point{};
    while (true) {    // only accepts point to be inside the sphere
        point = pickInCuboid<DIM>(randGenerator, randomReal, length);
        for (auto& coord : point) {
            coord -= 0.5;
        }
        if (geomTools::normeCarre<DIM>(point) < 0.5 * 0.5) {
            break;
        }
    }
    geomTools::renormalize<DIM>(point);
    return point;
}
} // namespace randomShooter
} // namespace sac_de_billes

#endif /* RANDOMSHOOTER_IXX_ */
