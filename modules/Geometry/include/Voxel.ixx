//! Copyright : see license.txt
//!
//! \brief

#pragma once


namespace merope {
namespace vox {
namespace auxi {
template<unsigned short DIM, class C>
Point<DIM> auxi::origin(const C& discPt,
    const Point<DIM>& dx) {
    Point <DIM> res;
    for (unsigned short i = 0; i < DIM; i++) {
        res[i] = discPt[i] * dx[i];
    }
    return res;
}

template<unsigned short DIM>
Point<DIM> auxi::center(const DiscPoint<DIM>& discPt,
    const Point<DIM>& dx) {
    Point <DIM> result;
    for (size_t i = 0; i < DIM; i++) {
        result[i] = (discPt[i] + 0.5) * dx[i];
    }
    return result;
}
}  // namespace  auxi
}  // namespace  vox
}  // namespace  merope
