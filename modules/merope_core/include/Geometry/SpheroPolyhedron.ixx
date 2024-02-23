//! Copyright : see license.txt
//!
//! \brief
//!
#ifndef ALGOPACKING_SRC_GEOMETRY_SPHEROPOLYHEDRON_IXX_
#define ALGOPACKING_SRC_GEOMETRY_SPHEROPOLYHEDRON_IXX_

#include "../../../AlgoPacking/src/Geometry/GeomTools_1.hxx"
#include "../Geometry/GeomTools.hxx"
#include "../MeropeNamespace.hxx"
#include "AuxiConvexPolyhedron.hxx"


namespace sac_de_billes {

template<unsigned short DIM>
inline SpheroPolyhedron<DIM>::SpheroPolyhedron(
    Identifier phase_, const vector<Point<DIM> >& vertices,
    const vector<vector<long> >& face_indices, double minkowskiRadius__) : SpheroPolyhedron<DIM>() {
    this->phase = phase_;
    this->minkowskiRadius_ = minkowskiRadius__;
    //
    vector<Point<DIM>> renormalized_vertices = vertices;
    //
    this->center = average<DIM>(renormalized_vertices);
    for (auto& v : renormalized_vertices) {
        v = v - this->center;
    }
    //
    this->innerPolyhedron = merope::facesFromVertices<DIM>(renormalized_vertices, face_indices);
    //
    this->outerPolyhedron = this->innerPolyhedron;
    for (auto& face : this->outerPolyhedron) {
        face.c() += minkowskiRadius_;
    }
    ////
    for (size_t i = 0; i < this->innerPolyhedron.size(); i++) {
        const auto& single_face_indice = face_indices[i];
        size_t halfSize = single_face_indice.size();
        auto translationVector = minkowskiRadius_ * this->innerPolyhedron[i].vec();
        //
        vector<Point<DIM>> local_list_vertices(2 * halfSize);
        for (size_t j = 0; j < halfSize; j++) {
            local_list_vertices[j] = renormalized_vertices[single_face_indice[j]];
            local_list_vertices[halfSize + j] = local_list_vertices[j] + translationVector;
        }
        //
        vector<vector<long>> local_face_indices(halfSize + 2); // face_indices local to the added polyhedron on the boundary
        for (long j = 0; j < static_cast<long>(halfSize); j++) {
            local_face_indices[0].push_back(halfSize - 1 - j); // inverse rotation
            local_face_indices[1].push_back(halfSize + j);
            local_face_indices[2 + j] = { j, j + 1, static_cast<long>(auxi_function::fast_modulo(halfSize + j + 1, 2 * halfSize)), static_cast<long>(halfSize + j) };
        }
        this->boundaryPolyhedrons.push_back(merope::facesFromVertices<DIM>(local_list_vertices, local_face_indices));
    }
    //
    if constexpr (DIM == 3) {
        auto list_of_edges = merope::edgesFromVertices<DIM>(renormalized_vertices, face_indices);
        for (const auto& edge : list_of_edges) {
            this->boundaryCylinders.push_back(Cylinder<DIM>(edge, minkowskiRadius_));
        }
    }
    //
    for (const auto& vertex : renormalized_vertices) {
        this->boundarySpheres_.push_back(Sphere<DIM>(vertex, minkowskiRadius_, 0));
    }
}

template<unsigned short DIM>
inline double SpheroPolyhedron<DIM>::distanceTo(
    const Point<DIM>& vector_from_center_to_point) const {
    auto tangentSpace_ = this->closestTangentSpace(vector_from_center_to_point);
    return  tangentSpace_.distanceTo(vector_from_center_to_point);
}

template<unsigned short DIM>
inline bool SpheroPolyhedron<DIM>::guaranteeOutside(
    const Point<DIM>& vector_from_center_to_point,
    const double& halfDiagonal) const {
    for (const auto& hf : this->outerPolyhedron) {
        if (hf.distanceTo(vector_from_center_to_point) > halfDiagonal) {
            return false;
        }
    }
    return true;
}

template<unsigned short DIM>
inline bool SpheroPolyhedron<DIM>::isInside(
    const Point<DIM>& point) const {
    auto vector_from_center_to_point = point - center;
    if (insidePolyhedralPart(vector_from_center_to_point)) return true;
    if constexpr (DIM == 3) {
        if (geomTools::isInside_Union<DIM>(this->boundaryCylinders, vector_from_center_to_point)) return true;
    }
    if (geomTools::isInside_Union<DIM>(this->boundarySpheres_, vector_from_center_to_point)) return true;
    return false;
}

template<unsigned short DIM>
inline vector<Point<DIM> > SpheroPolyhedron<DIM>::getVertices() const {
    vector<Point<DIM>> theCenters;
    std::transform(this->boundarySpheres().begin(), this->boundarySpheres().end(), std::back_inserter(theCenters), [](const auto& sph) {
        return sph.center;
        });
    return theCenters;
}

template<unsigned short DIM>
inline bool SpheroPolyhedron<DIM>::insidePolyhedralPart(
    const Point<DIM>& vector_from_center_to_point) const {
    if (geomTools::isInside_Intersection<DIM>(this->innerPolyhedron, vector_from_center_to_point)) return true;
    for (const auto& poly : this->boundaryPolyhedrons) {
        if (geomTools::isInside_Intersection<DIM>(poly, vector_from_center_to_point)) return true;
    }
    return false;
}

template<unsigned short DIM>
inline HalfSpace<DIM> SpheroPolyhedron<DIM>::closestTangentSpace(
    Point<DIM> vector_from_center_to_point) const {
    // strategy based on the following facts
    // 1° The closest point X of the spheroPolyhedron exists. Its tangent space is orthogonal to [X, vector_from_center_to_point]
    // 2° For each component (=Halsfpace, cylinders, spheres) find the candidates (=points X the normal of which is colinear to [X, vector_from_center_to_point])
    // 3° For each candidate point, test whether it is on the surface of the spheropolyhedron or inside it
    //      (It may happen for sphere, for example, where there are 2 candidates, and at least 1 is inside the spheropolyhedron)
    // 4° among all these candidates, find the closest one
    vector<array<Point<DIM>, 2>> vec_candidates_normal{}; // store {candidate point, normal}
    for (const auto& hf : this->outerPolyhedron) {
        vec_candidates_normal.push_back(array<Point<DIM>, 2>{geomTools::projection<DIM>(hf, vector_from_center_to_point), hf.vec()});
    }
    if constexpr (DIM == 3) {
        for (const auto& cyl : this->boundaryCylinders) {
            Point<DIM> projOnAxis = geomTools::projection<DIM>(cyl.axis, vector_from_center_to_point);
            Point<DIM> normal = vector_from_center_to_point - projOnAxis; geomTools::renormalize<DIM>(normal);
            vec_candidates_normal.push_back(array<Point<DIM>, 2>{projOnAxis + cyl.radius * normal, normal});
            vec_candidates_normal.push_back(array<Point<DIM>, 2>{projOnAxis - cyl.radius * normal, -normal});
        }
    }
    for (const auto& sph : this->boundarySpheres_) {
        Point<DIM> normal = vector_from_center_to_point - sph.center; geomTools::renormalize<DIM>(normal);
        vec_candidates_normal.push_back({ sph.center + sph.radius * normal, normal });
        vec_candidates_normal.push_back({ sph.center - sph.radius * normal, -normal });
    }
    //
    std::sort(vec_candidates_normal.begin(), vec_candidates_normal.end(), [&vector_from_center_to_point](const auto& pair1, const auto& pair2) {
        return geomTools::distanceCarre<DIM>(pair1[0], vector_from_center_to_point) < geomTools::distanceCarre<DIM>(pair2[0], vector_from_center_to_point);
        });
    //
    for (auto& pair : vec_candidates_normal) {
        if (geomTools::isPointOnBoundary<DIM>(*this, this->center + pair[0], pair[1], geomTools::EPSILON * this->minkowskiRadius_)) {
            return HalfSpace<DIM>(pair[1], pair[0]);
        }
    }
    // no candidate has been found! error
    cerr << __PRETTY_FUNCTION__ << endl;
    cerr << "Minkowski radius = " << this->minkowskiRadius_ << endl;
    throw runtime_error("No tangent space has been found. There is a bug or a 0 Minkowski radius.");
}

} // namespace merope

#endif /* ALGOPACKING_SRC_GEOMETRY_SPHEROPOLYHEDRON_IXX_ */
