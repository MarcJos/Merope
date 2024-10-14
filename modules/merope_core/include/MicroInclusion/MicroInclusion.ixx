//! Copyright : see license.txt
//!
//! \brief
//

#pragma once

#include "../MeropeNamespace.hxx"

#include "../Geometry/Trim.hxx"


namespace merope {
namespace smallShape {
template<unsigned short DIM, class SOLID>
Cuboid<DIM> createCuboid(const SOLID& solid) {
    static_assert(is_same_v<vector<Point<DIM>>, SOLID>
        or is_same_v<Sphere<DIM>, SOLID>
        or is_same_v<Ellipse<DIM>, SOLID>
        or is_same_v<SpheroPolyhedron<DIM>, SOLID>
        or is_same_v<Cylinder<3>, SOLID>
        );
    // set of points
    if constexpr (is_same_v<vector<Point<DIM>>, SOLID>) {
        Point <DIM> x_min;
        x_min.fill({ 0.5 * numeric_limits<double>::max() });
        Point <DIM> x_max;
        x_max.fill({ -0.5 * numeric_limits<double>::max() });
        for (const auto& p : solid) {
            for (size_t i = 0; i < DIM; i++) {
                x_min[i] = min(p[i], x_min[i]);
                x_max[i] = max(p[i], x_max[i]);
            }
        }
        return Cuboid <DIM>(x_min, x_max);
    }
    // sphere
    else if constexpr (is_same_v<Sphere<DIM>, SOLID>) {
        Point <DIM> x_min;
        Point <DIM> x_max;
        const Point<DIM>& center = solid.center;
        const double& R = solid.radius;
        for (size_t i = 0; i < DIM; i++) {
            x_min[i] = center[i] - R;
            x_max[i] = center[i] + R;
        }
        return Cuboid <DIM>(x_min, x_max);
    }
    // ellipse
    else if constexpr (is_same_v<Ellipse<DIM>, SOLID>) {
        auto alphas = solid.getAlphas();
        return createCuboid<DIM>(Sphere<DIM>(solid.center, *(std::min_element(alphas.begin(), alphas.end())), 0));
    }
    // spheroPolyhedron
    else if constexpr (is_same_v<SpheroPolyhedron<DIM>, SOLID>) {
        return createCuboid<DIM>(solid.getVertices(), solid.getMinkowskiRadius());
    }
    // Cylinder
    else if constexpr (is_same_v<Cylinder<3>, SOLID>) {
        Point<DIM> normal = solid.axis[1] - solid.axis[0];
        double  l = geomTools::renormalize<DIM>(normal);
        size_t i = 0;
        double val = 0.8;
        Point<DIM> vec = create_array<DIM>(0.);
        for (size_t j = 0; j < DIM; j++) {
            if (abs(normal[j]) < val) {
                i = j;
                val = abs(normal[j]);
            }
        }
        vec[i] = 1;
        Point<DIM> e1 = geomTools::prodVec<DIM>(vec, normal);
        geomTools::renormalize<DIM>(e1);
        Point<DIM> e2 = geomTools::prodVec<DIM>(e1, normal);
        geomTools::renormalize<DIM>(e2);
        // cylinder parametrized by (normal, e1, e2)
        Point <DIM> x_min;
        Point <DIM> x_max;
        const Point<DIM>& center = geomTools::get_center<DIM>(solid);
        for (size_t i_ = 0; i_ < DIM; i_++) {
            Point<DIM> vecteur_dir = create_array<DIM>(0.);
            vecteur_dir[i_] = 1;
            double ps1 = geomTools::prodScal<DIM>(vecteur_dir, e1);
            double ps2 = geomTools::prodScal<DIM>(vecteur_dir, e2);
            double norme_ps = sqrt(ps1 * ps1 + ps2 * ps2);
            x_min[i_] = center[i_] - 0.5 * l * abs(normal[i_]) - solid.radius * norme_ps;
            x_max[i_] = center[i_] + 0.5 * l * abs(normal[i_]) + solid.radius * norme_ps;
        }
        return Cuboid <DIM>(x_min, x_max);
    }
}

template<unsigned short DIM, class SOLID>
inline Cuboid<DIM> createCuboid(const SOLID& solid, double margin) {
    auto smallCuboid = smallShape::createCuboid<DIM>(solid);
    for (auto i = 0; i < DIM; i++) {
        smallCuboid.x_min[i] -= margin;
        smallCuboid.x_max[i] += margin;
    }
    return smallCuboid;
}

/// smallShape::ConvexPolyhedron<DIM>

//! Polyhedron factory
template<unsigned short DIM>
inline ConvexPolyhedronInc<DIM> smallShape::PolyhedronFactory<DIM>::fromVertices(
    Identifier ident, const vector<Point<DIM>>& vertices, const vector<vector<long>>& face_indices) {
    vector<Point<DIM>> renormalized_vertices = vertices;
    //
    Point<DIM> center = average<DIM>(renormalized_vertices);
    //
    for (auto& v : renormalized_vertices) {
        v = v - center;
    }
    auto faces = facesFromVertices<DIM>(renormalized_vertices, face_indices);
    Cuboid<DIM> cuboid = smallShape::createCuboid<DIM>(renormalized_vertices);
    cuboid.x_min += center;
    cuboid.x_max += center;
    return ConvexPolyhedronInc<DIM>(ident, cuboid, center, sac_de_billes::ConvexPolyhedron<DIM>(center, faces));

}

template<unsigned short DIM>
inline SpheroPolyhedronInc<DIM> merope::smallShape::SpheroPolyhedronFactory<DIM>::fromVertices(
    Identifier phase, const vector<Point<DIM> >& vertices,
    const vector<vector<long> >& face_indices, double minkowskiRadius) {
    vector<Point<DIM>> renormalized_vertices = vertices;
    //
    Point<DIM> center = average<DIM>(renormalized_vertices);
    //
    for (auto& v : renormalized_vertices) {
        v = v - center;
    }
    auto faces = facesFromVertices<DIM>(renormalized_vertices, face_indices);
    Cuboid<DIM> cuboid = smallShape::createCuboid<DIM>(renormalized_vertices, minkowskiRadius);
    cuboid.x_min += center;
    cuboid.x_max += center;
    return SpheroPolyhedronInc<DIM>(phase, cuboid, center, sac_de_billes::SpheroPolyhedron<DIM>(phase, vertices, face_indices, minkowskiRadius));
}

template<unsigned short DIM, class SOLID>
inline long  MicroInclusion_<DIM, SOLID>::whichLayer(const Point<DIM>& point) const {
    for (size_t layerIndex = 0; layerIndex < getNbOfLayers(); layerIndex++) {
        if (not isInside(point, layerIndex)) {
            return layerIndex - 1;
        }
    }
    return getNbOfLayers() - 1;
}

//! Rectangle

template<unsigned short DIM>
inline Rectangle<DIM>::Rectangle(Identifier ident, Point<DIM> xmin,
    Point<DIM> xmax) : ConvexPolyhedron<DIM>(ident, createCuboid<DIM>(xmin, xmax), 0.5 * (xmin + xmax), { ConvexPolyhedron<DIM>(Cuboid<DIM>(xmin, xmax)) }) {}

template<unsigned short DIM, class SOLID>
inline void MicroInclusion_<DIM, SOLID>::linearTransform(
    const Point<DIM>& linTransform) {
    cuboid.linearTransform(linTransform);
    linearTransform::point <DIM>(this->center, linTransform);
    for (auto& solid : this->solids) {
        linearTransform::proceed<DIM>(solid, linTransform);
    }
}

template<unsigned short DIM, class SOLID>
inline void MicroInclusion_<DIM, SOLID>::enlarge(double layerWidth) {
    solids[0] = geomTools::enlarge<DIM>(solids[0], layerWidth);
    cuboid.enlarge(layerWidth);
    for (size_t i = 1; i < layerIncrement.size(); i++) {
        layerIncrement[i] += layerWidth;
    }
}

template<unsigned short DIM, class SOLID>
inline void MicroInclusion_<DIM, SOLID>::pushLayer(double layerWidth,
    PhaseType layerPhase_) {
    layerIncrement.push_back(layerWidth);
    layerPhases.insert(layerPhases.end() - 1, layerPhase_);
    auto newSolid = geomTools::trim<DIM>(solids[solids.size() - 1], layerWidth);
    if constexpr (is_same_v<SOLID, Sphere<DIM>> or is_same_v<SOLID, Ellipse<DIM>>) {
        newSolid.phase = layerPhase_;
    }
    solids.push_back(newSolid);
}

template<unsigned short DIM, class SOLID>
inline bool MicroInclusion_<DIM, SOLID>::guaranteeInside(
    const Point<DIM>& vector_from_center_to_point,
    const double& halfDiagonal, const size_t& layerIndex) const {
    return distanceTo(vector_from_center_to_point, layerIndex) < -halfDiagonal;
}

template<unsigned short DIM, class SOLID>
inline double MicroInclusion_<DIM, SOLID>::distanceTo(const Point<DIM>& vector_from_center_to_point,
    size_t layerIndex) const {
    return geomTools::distanceTo<DIM>(getInnerInclusions()[layerIndex], vector_from_center_to_point);
}


template<unsigned short DIM, class SOLID>
inline bool MicroInclusion_<DIM, SOLID>::isInside(const Point<DIM>& point,
    const size_t& layerIndex) const {
    assert(layerIndex < this->getNbOfLayers());
    return geomTools::isInside<DIM>(getInnerInclusions()[layerIndex], point);
}

template<unsigned short DIM, class SOLID>
inline MicroInclusion_<DIM, SOLID>::MicroInclusion_(const SOLID& solid) :
    MicroInclusion_<DIM, SOLID>(solid.phase, createCuboid<DIM>(solid), geomTools::get_center<DIM>(solid), solid) {}

template<unsigned short DIM, class SOLID>
inline bool MicroInclusion_<DIM, SOLID>::guaranteeOutside(
    const Point<DIM>& vector_from_center_to_point,
    const double& halfDiagonal, const size_t& layerIndex) const {
    if constexpr (is_same_v<SOLID, ConvexPolyhedron<DIM>>) {
        for (size_t i = 0; i < DIM; i++) {
            if (vector_from_center_to_point[i] + this->center[i] < this->cuboid.x_min[i] - halfDiagonal
                or vector_from_center_to_point[i] + this->center[i] > this->cuboid.x_max[i] + halfDiagonal) {
                return true;
            }
        }
    }
    return distanceTo(vector_from_center_to_point, layerIndex) > halfDiagonal;
}

template<unsigned short DIM, class SOLID>
inline PhaseType& MicroInclusion_<DIM, SOLID>::getPhaseGraphical(size_t i) {
    return layerPhases[auxi_MicroInclusions::getIndexPhaseGraphical(i, layerPhases.size())];
}

template<unsigned short DIM, class SOLID>
inline const PhaseType& MicroInclusion_<DIM, SOLID>::getPhaseGraphical(size_t i) const {
    return layerPhases[auxi_MicroInclusions::getIndexPhaseGraphical(i, layerPhases.size())];
}


}  // namespace smallShape
}  // namespace merope


