//! Copyright : see license.txt
//!
//! \brief
//
#pragma once

#include "../../GenericMerope/StdHeaders.hxx"

#include "GeomTypes.hxx"


namespace merope {

template<unsigned short DIM>
class SpheroPolyhedron {
    //! given a Polyhedron P, is the set { x \in \R^d, dist(x,P) <= Minkowski radius}
    //! \warning all the components of the spheroPolyhedron are computed in local coordinates (i.e. wrt to a reference point called "center")
    //! \warning most function should be used with vector_from_center_to_point = given_point - center
    static_assert(DIM == 2 or DIM == 3);
public:
    //! constructor
    SpheroPolyhedron(Identifier phase, const vector<Point<DIM>>& vertices, const vector<vector<long>>& face_indices, double minkowskiRadius);
    //! compute the signed distance between a point and the surface of the spheroPolyhedron
    double distanceTo(const Point<DIM>& vector_from_center_to_point) const;
    //! guarantee that the square with center vector_from_center_to_point and halfDiagonal fixed is outside of the spheroPolyhedron
    bool guaranteeOutside(const Point<DIM>& vector_from_center_to_point, const double& halfDiagonal) const;
    //! return whether the point is inside the spheropolyhedron
    bool isInside(const Point<DIM>& point) const;
    //! \return the tangent space to the spheroPolyhedron clostest to the point vector_from_center_to_point
    //! \warning : in local coordinates
    HalfSpace<DIM> closestTangentSpace(Point<DIM> vector_from_center_to_point) const;

    //! getter
    const vector<Sphere<DIM>>& boundarySpheres() const { return this->boundarySpheres_; }
    //! getter
    double getMinkowskiRadius() const { return this->minkowskiRadius_; }
    //! get the underlying vertices
    vector<Point<DIM>> getVertices() const;
    //! getter
    const vector<HalfSpace<DIM> >& getInnerPolyhedron() const { return innerPolyhedron; }
    //! getter
    const vector<HalfSpace<DIM> >& getOuterPolyhedron() const { return outerPolyhedron; }
    //! getter
    const Point<DIM>& getCenter() const { return center; }
    //! getter
    PhaseType getPhase()const { return phase; };

private:
    //! Minkowski radius.
    double minkowskiRadius_;
    //! center
    Point<DIM> center;
    //! phase
    PhaseType phase;

    //! inner polyhedron (=defined by initial vertices)
    vector<HalfSpace<DIM>> innerPolyhedron;
    //! outer polyhedron (=inner polyedron, with face translated along their normals by a distance = minkowskiRadius)
    vector<HalfSpace<DIM>> outerPolyhedron;
    //! border polyhedrons
    vector<vector<HalfSpace<DIM>>> boundaryPolyhedrons;
    //! if dimension = 3 borderCylidners
    typename std::conditional_t<DIM == 3, vector<Cylinder<3>>, bool> boundaryCylinders;
    //! border balls
    vector<Sphere<DIM>> boundarySpheres_;

    //! default constructor
    SpheroPolyhedron() :minkowskiRadius_{ 0. }, center{}, phase{ 0 }, innerPolyhedron{}, outerPolyhedron{}, boundaryPolyhedrons{}, boundaryCylinders{}, boundarySpheres_{} {}
    //! \return whether inside the polyhedral part of the spheropolyhedron
    bool insidePolyhedralPart(const Point<DIM>& vector_from_center_to_point) const;
};

}  // namespace merope

#include "SpheroPolyhedron.ixx"


