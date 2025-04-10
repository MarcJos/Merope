//! Copyright : see license.txt
//!
//! \brief

#pragma once


namespace merope {
namespace voroInterface {

template<>
inline VoroInterface<3>::VoroInterface(array<double, 3> L, const vector<Sphere<3>>& centerTessels, array<bool, 3> periodicity) :
    InsideTorus<3>(L),
    PreparedVoroppContainer(centerTessels, L, periodicity), virtualLength{ L } {
    if (*min_element(L.begin(), L.end()) < 1e-3) {
        std::cerr << "WARNING" << endl;
        std::cerr << __PRETTY_FUNCTION__ << endl;
        std::cerr << "It is NOT advised to give cube lengths not of order 1." << endl;
        std::cerr << "Indeed, voro++ does not behave well with very small lengths." << endl;
        std::cerr << "WARNING" << endl;
    }
}

template<>
inline VoroInterface<2>::VoroInterface(array<double, 2> L,
    const vector<Sphere<2>>& centerTessels, array<bool, 2> periodicity) :
    InsideTorus<2>(L),
    PreparedVoroppContainer(voroInterface_aux::extendDimension(centerTessels, voroInterface_aux::virtualLength<2>(L)), voroInterface_aux::virtualLength<2>(L), voroInterface_aux::extendPeriodicity(periodicity))
    , virtualLength{ voroInterface_aux::virtualLength<2>(L) } {
}


template<unsigned short DIM>
int VoroInterface<DIM>::findTessel(const Point<DIM>& pt) {
    if constexpr (DIM == 3) {
        return voroInterface_aux::findTessel(pt, voropp_container);
    } else if constexpr (DIM == 2) {
        return voroInterface_aux::findTessel(voroInterface_aux::extendDimension(pt, virtualLength), voropp_container);
    } else {
        Merope_assert(false, "Only in dimensions 2 and 3");
    }
}

template<unsigned short DIM>
void VoroInterface<DIM>::drawGnuPlot(string fileName) {
    voropp_container->draw_cells_gnuplot(fileName.c_str());
}

template<unsigned short DIM>
void VoroInterface<DIM>::drawCellsPov(string fileName) {
    voropp_container->draw_cells_pov(fileName.c_str());
}

template<unsigned short DIM>
void VoroInterface<DIM>::printCustom(string format, string fileName) {
    voropp_container->print_custom(format.c_str(), fileName.c_str());
}

template<unsigned short DIM>
vector<smallShape::ConvexPolyhedronInc<DIM>> VoroInterface<DIM>::getMicroInclusions() {
    vector<smallShape::ConvexPolyhedronInc<DIM>> polyhedrons{ };
    polyhedrons.reserve(voropp_container->total_particles());
    //
    auto addInclusion_ = [&](auto pp, auto index, auto c_pt) {
        addInclusion(pp, index, c_pt, polyhedrons);
        };
    loop_on_voroppcontainer(voropp_container, addInclusion_);
    return polyhedrons;
}

template<unsigned short DIM>
std::pair<std::map<Identifier, vector<Identifier>>, std::map<Identifier, HalfSpace<3>>> VoroInterface<DIM>::computeSolids() {
    return auxi::computeSolids(voroInterface_aux::virtualLength<DIM>(this->getL()), this->getSingleCells(), this->periodicity);
}

template<unsigned short DIM>
std::map<Identifier, Point<3>> VoroInterface<DIM>::getCellCenters() {
    if constexpr (DIM == 2) {
        Merope_error_not_done();
    }
    std::map<Identifier, Point<3>> centers;
    vector<merope::voroInterface::SingleCell> polys = this->getSingleCells();

    for (Identifier id_solid = 0; id_solid < polys.size(); id_solid++) {
        centers.insert({ polys[id_solid].identifier, polys[id_solid].center });
    }
    return centers;
}

template<unsigned short DIM>
inline vector<SingleCell> VoroInterface<DIM>::getSingleCells() {
    if constexpr (DIM == 2) {
        Merope_error_not_done();
    }
    vector<SingleCell> polyhedrons{ };
    polyhedrons.reserve(voropp_container->total_particles());
    //
    auto addCell_ = [&](auto pp, auto index, auto c_pt) {
        addCell(pp, index, c_pt, polyhedrons);
        };
    loop_on_voroppcontainer(voropp_container, addCell_);
    return polyhedrons;
}


template<unsigned short DIM>
void VoroInterface<DIM>::addWallCylinder(double xc_, double yc_, double zc_, double xa_, double ya_, double za_, double rc_) {
    PreparedVoroppContainer::addWallCylinder(xc_, yc_, zc_, xa_, ya_, za_, rc_);
}

template<unsigned short DIM>
void VoroInterface<DIM>::addInclusion(double* pp, int index,
    voro::voronoicell_neighbor* c,
    vector<smallShape::ConvexPolyhedronInc<DIM>>& polyhedrons) {
    vector<HalfSpace<DIM>> faces = voroInterface_aux::getFaces<DIM>(c);  // the faces
    Point<DIM> center;  // center of the polyhedron
    if constexpr (DIM == 3) {
        center = { pp[0], pp[1], pp[2] };
    } else  if constexpr (DIM == 2) {
        center = { pp[0], pp[1] };
    }
    Cuboid<DIM> cuboid = voroInterface_aux::getCuboid<DIM>(c);  // the surrounding cuboid
    //
    for (size_t i = 0; i < DIM; i++) {
        cuboid.x_min[i] += center[i];
        cuboid.x_max[i] += center[i];
    }
    //
    polyhedrons.push_back(
        smallShape::ConvexPolyhedronInc<DIM>(index, cuboid, center, sac_de_billes::ConvexPolyhedron<DIM>(center, faces)));
}

template<unsigned short DIM>
void VoroInterface<DIM>::addCell(double* pp, int index,
    voro::voronoicell_neighbor* c, vector<SingleCell>& polyhedrons) {
    Point<3> center;  // center of the polyhedron
    if constexpr (DIM == 3) {
        center = { pp[0], pp[1], pp[2] };
    } else  if constexpr (DIM == 2) {
        center = { pp[0], pp[1], 0 };
    }
    polyhedrons.emplace_back(SingleCell(index, center, c));
}



// voroInterface_aux
template<unsigned short DIM>
Point<3> voroInterface_aux::virtualLength(Point<DIM> L) {
    if constexpr (DIM == 2) {
        return Point<3> {L[0], L[1],
            VIRTUAL_LENGTH_MULTIPLIER* min(L[0], L[1])};
    } else if constexpr (DIM == 3) {
        return L;
    }
}


template<unsigned short DIM>
Cuboid<DIM> voroInterface_aux::getCuboid(voro::voronoicell_neighbor* cell) {
    if (not cell) {
        return Cuboid<DIM>(create_array<DIM>(0.), create_array<DIM>(0.));
    }
    //
    static_assert(DIM == 2 or DIM == 3);
    Cuboid<3> cuboid3D = smallShape::createCuboid<3>(getRenormalizedVertices(cell));
    if constexpr (DIM == 3) {
        return cuboid3D;
    } else if constexpr (DIM == 2) {
        Point<2> x_min = geomTools::linearProjection<2, 3>(cuboid3D.x_min);
        Point<2> x_max = geomTools::linearProjection<2, 3>(cuboid3D.x_max);
        return Cuboid<2>(x_min, x_max);
    }
}

template<unsigned short DIM>
inline vector<HalfSpace<DIM> > voroInterface_aux::getFaces(
    voro::voronoicell_neighbor* cell) {
    static_assert(DIM == 2 or DIM == 3);
    if constexpr (DIM == 2)        return getFaces2D(cell);
    else if constexpr (DIM == 3)   return getFaces3D(cell);
}

template<unsigned short DIM>
inline std::map<Identifier, vector<Identifier>> getClosestNeighbors(
    array<double, DIM> L, const vector<Sphere<DIM> >& centerTessels) {
    voroInterface::VoroInterface<DIM> voroInt(L, centerTessels);
    auto cells = voroInt.getSingleCells();
    std::map<Identifier, vector<Identifier>> result{};
    for (const auto& cell : cells) {
        result[cell.identifier] = cell.neighbors;
    }
    return result;
}

template<class Func>
void loop_on_voroppcontainer(voro::container_poly* voropp_container, Func my_function) {
    voro::c_loop_all vl(*voropp_container);
    // Heavily inspired from voro::print_custom
    int ijk, q;
    double* pp;
    voro::voronoicell_neighbor c;
    if (vl.start()) do {
        ijk = vl.ijk;
        q = vl.q;
        int index = voropp_container->id[ijk][q];  // index of the polyhedron
        pp = voropp_container->p[ijk] + voropp_container->ps * q;
        if (voropp_container->compute_cell(c, vl)) {
            my_function(pp, index, &c);
        } else {
            my_function(pp, index, nullptr);
        }
    } while (vl.inc());
}

namespace auxi {
template<bool Ignore_Negative_Neighbor>
pair<Identifier, vector<size_t>> correspondingFaces(std::function<const merope::voroInterface::SingleCell* (Identifier)> getCell,
    Identifier id_cell, size_t id_face) {
    const auto& singleCell_0 = getCell(id_cell);
    const auto& normal_0 = singleCell_0->faceNormal[id_face];
    //
    Identifier id_singleCell_1 = singleCell_0->neighbors[id_face];
    // verifications
    if (id_singleCell_1 < 0) {
        if constexpr (Ignore_Negative_Neighbor) {
            return { 0,{} };
        } else {
            Merope_assert(false, "Voro++ is undefined with negative neighbors (at least, does not correspond to the users' manual)");
        }
    }
    //
    const voroInterface::SingleCell* singleCell_1 = getCell(id_singleCell_1);
    //
    if (singleCell_1) {
        vector<size_t> candidatesFaces = {};
        double tol = 1e-16; // magical constant!
        for (size_t face_1_id_ = 0; face_1_id_ < singleCell_1->faceVertices.size(); face_1_id_++) {
            if (geomTools::normeCarre<3>(normal_0 + singleCell_1->faceNormal[face_1_id_]) < tol) {
                if (singleCell_1->neighbors[face_1_id_] == id_cell) {
                    candidatesFaces.push_back(face_1_id_);
                } else {
                    if constexpr (Ignore_Negative_Neighbor) {
                        Merope_assert(singleCell_1->neighbors[face_1_id_] + id_cell == 0, "Two faces have close normals, but do not correspond to each others");
                    } else {
                        std::cerr << singleCell_1->neighbors[face_1_id_] << " : " << id_cell << endl;
                        Merope_assert(false, "Two faces have close normals, but do not correspond to each others");
                    }
                }
            }
        }
        return pair<Identifier, vector<size_t>>(id_singleCell_1, candidatesFaces);
    } else {
        return { 0,{} };
    }
}

}  // namespace  auxi
}  // namespace voroInterface
}  // namespace merope



