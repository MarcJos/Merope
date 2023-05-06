//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef POINT_IXX_
#define POINT_IXX_

namespace sac_de_billes {

/*
template<unsigned short DIM>
inline Point<DIM>::Point(const RenormPoint<DIM>& renormPoint): array<double, DIM>(renormPoint.getPoint()) {}
*/

template<unsigned short DIM>
static ostream& print(ostream& out, const Point<DIM>& p) {
    out << "Point: ";
    out << p[0];
    for (unsigned short i = 1; i < DIM; i++) {
        out << ",  " << p[i];
    }
    out << "\n";
    return out;
}

template<unsigned short DIM>
ostream& operator<<(ostream& out, const Point<DIM>& p) {
    return print(&out, &p);
}

template<unsigned short DIM, class C>
constexpr array<C, DIM> create_array(C x) {
    // fixme
    static_assert(DIM == 1 or DIM == 2 or DIM == 3);
    if constexpr (DIM == 1) {
        return array<C, DIM> {x};
    }
    else if constexpr (DIM == 2) {
        return array<C, DIM>{x, x};
    }
    else if constexpr (DIM == 3) {
        return array<C, DIM>{x, x, x};
    }
}

template<unsigned short DIM>
inline Point<DIM> average(const vector<Point<DIM>>& vPoints) {
    Point<DIM> result = create_array<DIM>(0.);
    for (const auto& pt : vPoints) {
        result += pt;
    }
    result *= 1. / vPoints.size();
    return result;
}

} // namespace sac_de_billes

#endif /* POINT_IXX_ */
