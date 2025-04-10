//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

namespace sac_de_billes {

namespace randomShooter {

template<unsigned short DIM>
inline Point<DIM> pickInCuboid(mt19937& randGenerator,
    const Point<DIM>& length) {
    static uniform_real_distribution<> randomReal(0., 1.);
    //
    Point <DIM> point{};
    for (size_t i = 0; i < DIM; i++) {
        point[i] = randomReal(randGenerator) * length[i];
    }
    return point;
}

template<unsigned short DIM>
inline Point<DIM> pickOnSphere(mt19937& randGenerator) {
    //
    Point<DIM> length = create_array<DIM>(1.);
    Point<DIM> point{};
    while (true) {    // only accepts point to be inside the sphere
        point = pickInCuboid<DIM>(randGenerator, length);
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
}  // namespace randomShooter
}  // namespace sac_de_billes


