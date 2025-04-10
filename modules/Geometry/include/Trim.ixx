//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


template<unsigned short DIM, class SOLID>
inline SOLID merope::geomTools::trim(const SOLID& solid,
    double layerWidth) {
    static_assert(is_same_v<SOLID, Sphere<DIM>> or is_same_v<SOLID, ConvexPolyhedron<DIM>>
        or is_same_v<SOLID, Cylinder<3>> or is_same_v<SOLID, SpheroPolyhedron<DIM>>
        or is_same_v<SOLID, Ellipse<DIM>>);
    // case 1 : sphere
    if constexpr (is_same_v<SOLID, Sphere<DIM>>) {
        return Sphere<DIM>(solid.center, solid.radius - layerWidth, solid.phase);
    }
    // case 2 : convex polyhedron
    else if constexpr (is_same_v<SOLID, ConvexPolyhedron<DIM>>) {
        vector<HalfSpace<DIM>> faces = solid.faces;
        for (auto& hs : faces) {
            hs.c() -= layerWidth;
        }
        return ConvexPolyhedron<DIM>(solid.center, faces);
    }
    // case 3 : cylinder
    else if constexpr (is_same_v<SOLID, Cylinder<3>>) {
        auto res = solid;
        res.radius -= layerWidth;
        auto vec = get<1>(res.axis) - get<0>(res.axis);
        geomTools::renormalize<DIM>(vec);
        get<0>(res.axis) += layerWidth * vec;
        get<1>(res.axis) -= layerWidth * vec;
        return res;
    }
    // case 4 : spheroPolyhedron
    else if constexpr (is_same_v<SOLID, SpheroPolyhedron<DIM>>) {
        double minkowski_radius = solid.getMinkowskiRadius();
        if (layerWidth < minkowski_radius) {
            throw runtime_error("Correct internal spheropolyhedron should be computed");
        } else {
            throw runtime_error("");
            /* Correct answer : but incorrect type
            vector<HalfSpace<DIM>> faces = solid.getInnerPolyhedron();
            for (auto& hs : faces) {
                hs.c() -= layerWidth - minkowski_radius;
            }
            return ConvexPolyhedron<DIM>(solid.center, faces);
            */
        }
    }
    // case 5 : ellipse
    else if constexpr (is_same_v<SOLID, Ellipse<DIM>>) {
        Merope_error_not_done();
    }
    //
    else {
        Merope_static_error(SOLID, "Impossible");
    }
}


template<unsigned short DIM, class SOLID>
SOLID merope::geomTools::enlarge(const SOLID& solid, double layerWidth) {
    static_assert(is_same_v<SOLID, Sphere<DIM>> or is_same_v<SOLID, ConvexPolyhedron<DIM>>
        or is_same_v<SOLID, Cylinder<3>> or is_same_v<SOLID, SpheroPolyhedron<DIM>>
        or is_same_v<SOLID, Ellipse<DIM>>);
    // case 1 : sphere
    if constexpr (is_same_v<SOLID, Sphere<DIM>>) {
        return Sphere<DIM>(solid.center, solid.radius + layerWidth, solid.phase);
    }
    // case 2 : convex polyhedron
    else if constexpr (is_same_v<SOLID, ConvexPolyhedron<DIM>>) {
        vector<HalfSpace<DIM>> faces = solid.faces;
        for (auto& hs : faces) {
            hs.c() += layerWidth;
        }
        return ConvexPolyhedron<DIM>(solid.center, faces);
    }
    // case 3 : cylinder
    else if constexpr (is_same_v<SOLID, Cylinder<3>>) {
        auto res = solid;
        res.radius += layerWidth;
        auto vec = get<1>(res.axis) - get<0>(res.axis);
        geomTools::renormalize<DIM>(vec);
        get<0>(res.axis) -= layerWidth * vec;
        get<1>(res.axis) += layerWidth * vec;
        return res;
    }
    // case 4 : spheroPolyhedron
    else if constexpr (is_same_v<SOLID, SpheroPolyhedron<DIM>>) {
        Merope_error_not_done();
    }
    // case 5 : ellipse
    else if constexpr (is_same_v<SOLID, Ellipse<DIM>>) {
        Merope_error_not_done();
    }
    //
    else {
        Merope_static_error(SOLID, "Impossible");
    }
}


