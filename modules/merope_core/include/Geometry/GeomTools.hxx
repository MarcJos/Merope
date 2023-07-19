//! Copyright : see license.txt
//!
//! \brief 
//!

#ifndef GEOMTOOLS_HXX_
#define GEOMTOOLS_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../../../AlgoPacking/src/AmbiantSpace.hxx"


namespace sac_de_billes {
namespace geomTools {

enum class Intersection_LineConvex {
        //! 4 possibilities of intersection between a line (x_0, y_0, z) and a convexShape
        Empty, // no intersection
        Plus,  // z > z_0
        Minus, // z < z_0
        All,    // all z
        Segment // z_0 < z < z_1
};

//! for describing 3 possibilities
enum class TypeFace {
        Left = -1, Right = 1, None = 0
};

//! \return the projection of a point on the plane delimitating the Half-space
//! \param point : point to be projected
//! \param halfspace : describes the plan on which we project
template<unsigned short DIM>
Point<DIM> projection(const HalfSpace<DIM>& halfspace, const Point<DIM>& point);
//! \return the closest Half-space to a point
//! \param point : point to be projected
//! \param halfspaces : describes the polyhedron on which we project
template<unsigned short DIM>
const HalfSpace<DIM>& closestHalfplane(const vector<HalfSpace<DIM>>& halfplaneList, Point<DIM> x);


//! \return the equation of a given face of a cuboid
//! \param cuboid
//! \param direction : 0, 1, 2 for directions x, y, z
//! \param typeFace : left or right face
//! \warning  : not tested
template<unsigned short DIM>
HalfSpace<DIM> faceEquation(const AmbiantSpace::Cube<DIM>& cuboid, int direction, TypeFace typeFace);

//! \return the point of intersection between a line and a plane (assume the line is not parallel to the plane)
//! \param halfspace : plane
//! \param A, B : extremities of the line
//! \warning  : not tested
template<unsigned short DIM>
Point<DIM> intersectLinePlane(const HalfSpace<DIM>& halfSpace, const Point<DIM>& A, const Point<DIM>& B);

//! computes the intersection with a line (x1, x2, z) or (x1, z), where z is in |R and the half-space
//! \return the type of intersection, put the result in z, with the sign of the halfline
//! fixme
template<unsigned short DIM>
Intersection_LineConvex computeIntersection(const HalfSpace<DIM>& halfSpace, const array<double, DIM - 1>& x1x2, double& z,
        double distance = 0);
//! computes the intersection with a line (x1, x2, z) or (x1, z) and the half-space
//! \return the type of intersection, else put the result in z, with the sign of the halfline
//! \param SIGN : upper or lower intersection
template<unsigned short DIM, unsigned short SIGN>
Intersection_LineConvex computeIntersection(const Sphere<DIM>& sphere, const array<double, DIM - 1>& x1x2, double& z);

//! \return the volume of the intersection of a unit cube with a half-space
//! \see The Volume of Simplices Clipped by a Half Space, CHO and CHO, Applied Mathematics Letters, 2000
template<unsigned short DIM>
double fracVolIntersection(const HalfSpace<DIM>& hf);
//! additional parameter here : the cuboid is supposed to be of form [0,cubeLenght[0]] x [0,cubeLenght[1]] ...
//! the result is the volume fraction of cuboid inside the half-space, wrt to the full cuboid
template<unsigned short DIM>
double fracVolIntersection(const Point<DIM>& cubeLength, HalfSpace<DIM> hf);
//! for a unique face
template<unsigned short DIM>
double fracVolIntersection(HalfSpace<DIM> face,
        const Point<DIM>& centerPoly_to_origVoxel,
        const Point<DIM>& cubeLength);
//! for a list of faces
template<unsigned short DIM>
double fracVolIntersection(const vector<HalfSpace<DIM>>& faces,
        const Point<DIM>& centerPoly_to_origVoxel,
        const Point<DIM>& cubeLength);
//! for a spheroPolyhedron
template<unsigned short DIM>
double fracVolIntersection(const SpheroPolyhedron<DIM>& sphP,
        const Point<DIM>& centerPoly_to_origVoxel,
        const Point<DIM>& cubeLength);

namespace fracVolIntersect {
//! \brief auxiliary functions for geomTools::fracVolIntersection
template<unsigned short DIM>
double core(const HalfSpace<DIM>& hf);
template<unsigned short DIM>
double auxi(const HalfSpace<DIM>& hf);
}

} // namespace geomTools

namespace auxi_Corner_of_Cubes {
template<unsigned short DIM>
constexpr size_t NumberOfCorners() {
        if constexpr (DIM == 2)
                return 4;
        else if constexpr (DIM == 3)
                return 8;
        else throw logic_error(__PRETTY_FUNCTION__);
}

//! stores the result of TABCORNER
constexpr array<array<double, 3>, 8> TABCORNER3D{ array<double, 3> { 0, 0, 0 },
        array<double, 3> { 1, 1, 1 }, array<double, 3> { 0, 0, 1 }, array<
                double, 3> { 1, 1, 0 }, array<double, 3> { 0, 1, 0 }, array<
                double, 3> { 1, 0, 1 }, array<double, 3> { 0, 1, 1 }, array<
                double, 3> { 1, 0, 0 } };
constexpr array<array<double, 2>, 4> TABCORNER2D{ array<double, 2> { 0, 0 },
        array<double, 2> { 1, 1 }, array<double, 2> { 0, 1 }, array<double, 2> {
                1, 0 } };
//! stores the result of Indices_TABCORNER
constexpr array<short, 8> Indices_TABCORNER3D{ -1, 1, 1, -1, 1, -1, -1, 1 };
constexpr array<short, 4> Indices_TABCORNER2D{ 1, 1, -1, -1 };
} // namespace auxi_Corner_of_Cubes

namespace Corners_of_Cubes {
//! Class to manage corners of cubes
//! opposite corners in order to reach get more easily out of the sphere
template<unsigned short DIM>
constexpr const array<array<double, DIM>,
        auxi_Corner_of_Cubes::NumberOfCorners<DIM>()>& TabCorner();
//! indicates for each TAB, the number of 0 per corner
template<unsigned short DIM>
constexpr const array<short, auxi_Corner_of_Cubes::NumberOfCorners<DIM>()>& Indices_TabCorner();

} // namespace Corner_of_Cubes
} // namespace sac_de_billes


#include "../Geometry/GeomTools.ixx"

#endif /* GEOMTOOLS_HXX_ */
