//! Copyright : see license.txt
//!
//! \brief 
//!

#ifndef GEOMTOOLS_IXX_
#define GEOMTOOLS_IXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"


namespace sac_de_billes {

// free functions
template<unsigned short DIM>
inline Point<DIM> geomTools::projection(const HalfSpace<DIM>& halfspace, const Point<DIM>& point) {
    return point - (geomTools::prodScal<DIM>(point, halfspace.vec()) - halfspace.c()) * halfspace.vec();
}



template<unsigned short DIM>
inline const HalfSpace<DIM>& geomTools::closestHalfplane(
    const vector<HalfSpace<DIM> >& polyhedron, Point<DIM> point) {
    return *(std::min_element(polyhedron.begin(), polyhedron.end(), [&point](const auto& hf1, const auto& hf2) {
        return hf1.distanceTo(point) < hf2.distanceTo(point);
        }));
}

template<unsigned short DIM>
inline HalfSpace<DIM> geomTools::faceEquation(
    const AmbiantSpace::Cube<DIM>& cuboid, int direction,
    TypeFace typeFace) {
    if (typeFace == TypeFace::None) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Unexpected");
    }
    double c = (typeFace == TypeFace::Left) ? 0 : cuboid.L[direction];
    Point<DIM> xi = create_array<DIM>(0.);
    xi[direction] = (typeFace == TypeFace::Left) ? -1 : 1;
    return HalfSpace<DIM>(xi, c);
}

template<unsigned short DIM>
inline Point<DIM> geomTools::intersectLinePlane(const HalfSpace<DIM>& halfSpace,
    const Point<DIM>& A, const Point<DIM>& B) {
    Point<DIM> AB = B - A;
    double t = (halfSpace.c() - geomTools::prodScal<DIM>(A, halfSpace.vec())) / (geomTools::prodScal<DIM>(AB, halfSpace.vec()));
    return A + (t * AB);
}

template<unsigned short DIM>
inline geomTools::Intersection_LineConvex geomTools::computeIntersection(const HalfSpace<DIM>& halfSpace, const array<double, DIM - 1>& x1x2, double& z, double distance) {
    double sumOther = halfSpace.c() + distance;
    for (size_t i = 0; i < DIM - 1; i++) {
        sumOther -= halfSpace.vec()[i] * x1x2[i];
    }
    //! fixme : magical constant
    if (abs(halfSpace.vec()[DIM - 1]) < 1.e-10) { // consider that last coordinate of vec is 0
        z = numeric_limits<double>::max();;
        if (sumOther > 0) {
            return Intersection_LineConvex::All;
        }
        else {
            return Intersection_LineConvex::Empty;
        }
    }
    else {
        z = sumOther / halfSpace.vec()[DIM - 1];
        if (halfSpace.vec()[DIM - 1] > 0) {
            return Intersection_LineConvex::Minus;
        }
        else {
            return Intersection_LineConvex::Plus;
        }
    }
}


template<unsigned short DIM, unsigned short SIGN>
inline geomTools::Intersection_LineConvex geomTools::computeIntersection(
    const Sphere<DIM>& sphere, const array<double, DIM - 1>& x1x2,
    double& z) {
    double heightWrtCenter = auxi_function::puissance<2>(sphere.radius) - geomTools::normeCarre<DIM - 1>(x1x2);
    // first case : emplty intersection
    if (heightWrtCenter < 0) {
        return geomTools::Intersection_LineConvex::Empty;
    }
    // second case : a segment
    heightWrtCenter = sqrt(heightWrtCenter);
    z = sphere.center[DIM - 1] + SIGN * heightWrtCenter;
    return geomTools::Intersection_LineConvex::Segment;
}

///fracVolIntersection

template<unsigned short DIM>
inline double geomTools::fracVolIntersection(
    const vector<HalfSpace<DIM>>& faces,
    const Point<DIM>& centerPoly_to_origVoxel,
    const Point<DIM>& cubeLength) {
    double volFrac = 1.; // volume fraction
    Point <DIM> linTransform; // inverse of cubeLength
    for (size_t i = 0; i < DIM; i++) {
        linTransform[i] = 1. / cubeLength[i];
    }
    for (auto face : faces) {
        face.c() -= geomTools::prodScal<DIM>(centerPoly_to_origVoxel, face.vec()); // the center of the voxel becomes the new origin
        linearTransform::proceed <DIM>(face, linTransform);
        volFrac *= geomTools::fracVolIntersection<DIM>(face);
        if (volFrac < geomTools::EPSILON) {
            return 0.;
        }
    }
    return volFrac;
}

template<unsigned short DIM>
inline double geomTools::fracVolIntersection(HalfSpace<DIM> face,
    const Point<DIM>& centerPoly_to_origVoxel,
    const Point<DIM>& cubeLength) {
    face.c() -= geomTools::prodScal<DIM>(centerPoly_to_origVoxel, face.vec()); // the center of the voxel becomes the new origin
    return geomTools::fracVolIntersection<DIM>(cubeLength, face);
}

template<unsigned short DIM>
inline double sac_de_billes::geomTools::fracVolIntersection(
    const SpheroPolyhedron<DIM>& sphP,
    const Point<DIM>& centerPoly_to_origVoxel,
    const Point<DIM>& cubeLength) {
    return fracVolIntersection<DIM>(sphP.closestTangentSpace(centerPoly_to_origVoxel + 0.5 * cubeLength), centerPoly_to_origVoxel, cubeLength);
}

template<unsigned short DIM>
inline double geomTools::fracVolIntersection(const Point<DIM>& cubeLength,
    HalfSpace<DIM> halfspace) {
    Point <DIM> linTransform; // inverse of cubeLength
    for (size_t i = 0; i < DIM; i++) {
        linTransform[i] = 1. / cubeLength[i];
    }
    linearTransform::proceed <DIM>(halfspace, linTransform);
    return fracVolIntersection <DIM>(halfspace);
}

template<unsigned short DIM>
inline double geomTools::fracVolIntersection(const HalfSpace<DIM>& hspace) {
    // tests if outerNormal close to 0
    const Point<DIM>& outerNormal = hspace.vec();
    // test if obviously outside or inside
    constexpr double sqrtDIM = sqrt(DIM);
    if (hspace.c() > sqrtDIM) {
        return 1;
    }
    else if (hspace.c() < -sqrtDIM) {
        return 0;
    }
    // proceed with numerical computation
    double norme2 = normeCarre <DIM>(outerNormal);
    if constexpr (DIM == 3) {
        if (abs(auxi_function::puissance < 2 >(outerNormal[0]))
            < EPSILON * norme2) {
            return fracVolIntersection < 2 >(HalfSpace<2> { Point<2> {
                outerNormal[1], outerNormal[2] }, hspace.c() });
        }
        else if (abs(auxi_function::puissance < 2 >(outerNormal[1]))
            < EPSILON * norme2) {
            return fracVolIntersection < 2 >(HalfSpace<2> { Point<2> {
                outerNormal[0], outerNormal[2] }, hspace.c() });
        }
        else if (abs(auxi_function::puissance < 2 >(outerNormal[2]))
            < EPSILON * norme2) {
            return fracVolIntersection < 2 >(HalfSpace<2> { Point<2> {
                outerNormal[0], outerNormal[1] }, hspace.c() });
        }
    }
    else if constexpr (DIM == 2) {
        if (abs(auxi_function::puissance < 2 >(outerNormal[0]))
            < EPSILON * norme2) {
            return fracVolIntersection < 1 >(HalfSpace<1> { Point<1> {
                outerNormal[1] }, hspace.c() });
        }
        else if (abs(auxi_function::puissance < 2 >(outerNormal[1]))
            < EPSILON * norme2) {
            return fracVolIntersection < 1 >(HalfSpace<1> { Point<1> {
                outerNormal[0] }, hspace.c() });
        }
    }
    // computes the intersection
    double volume = fracVolIntersect::core <DIM>(hspace);
    return volume;
}

// geomTools::fracVolIntersect::

template<unsigned short DIM>
inline double geomTools::fracVolIntersect::core(const HalfSpace<DIM>& hspace) {
    const Point<DIM>& outerNormal = hspace.vec();
    // case 1
    if constexpr (DIM == 1) {
        if (outerNormal[0] > 0) {
            return max(min(1., hspace.c() / outerNormal[0]), 0.);
        }
        else {
            return max(min(1., 1. - hspace.c() / outerNormal[0]), 0.);
        }
    }
    else {
        // case 2 or 3
        double prefactor;
        if constexpr (DIM == 2) {
            prefactor = 1. / (2 * outerNormal[0] * outerNormal[1]);
        }
        else if constexpr (DIM == 3) {
            prefactor = -1.
                / (6 * outerNormal[0] * outerNormal[1] * outerNormal[2]);
        }
        return prefactor * fracVolIntersect::auxi <DIM>(hspace);
    }
}

template<unsigned short DIM>
inline double geomTools::fracVolIntersect::auxi(const HalfSpace<DIM>& hspace) {
    double volume = 0;
    for (size_t i = 0; i < Corners_of_Cubes::TabCorner<DIM>().size(); i++) {
        double temp = hspace.c() - prodScal < DIM
        >(hspace.vec(), Corners_of_Cubes::TabCorner<DIM>()[i]);
        temp = (temp > 0) ? temp : 0; //positive part
        volume += Corners_of_Cubes::Indices_TabCorner<DIM>()[i]
            * auxi_function::puissance <DIM>(temp);
    }
    return volume;
}

// Corners of Cubes

template<unsigned short DIM>
constexpr const array<array<double, DIM>,
    auxi_Corner_of_Cubes::NumberOfCorners<DIM>()>& Corners_of_Cubes::TabCorner() {
    if constexpr (DIM == 2)
        return auxi_Corner_of_Cubes::TABCORNER2D;
    else if constexpr (DIM == 3)
        return auxi_Corner_of_Cubes::TABCORNER3D;
    else throw logic_error(__PRETTY_FUNCTION__);
}

template<unsigned short DIM>
constexpr const array<short, auxi_Corner_of_Cubes::NumberOfCorners<DIM>()>& Corners_of_Cubes::Indices_TabCorner() {
    if constexpr (DIM == 2)
        return auxi_Corner_of_Cubes::Indices_TABCORNER2D;
    else if constexpr (DIM == 3)
        return auxi_Corner_of_Cubes::Indices_TABCORNER3D;
    else throw logic_error(__PRETTY_FUNCTION__);
}

}//namespace sac_de_billes


#endif /* GEOMTOOLS_IXX_ */
