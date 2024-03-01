//! Copyright : see license.txt
//!
//! \brief

#ifndef VOROINTERFACE_IXX_
#define VOROINTERFACE_IXX_


#include "../MeropeNamespace.hxx"


namespace merope {
namespace voroInterface {

template<unsigned short DIM>
int VoroInterface<DIM>::findTessel(const Point<DIM>& pt) {
    if constexpr (DIM == 3) {
        return voroInterface_aux::findTessel(pt, voropp_container);
    } else if constexpr (DIM == 2) {
        return voroInterface_aux::findTessel(voroInterface_aux::extendDimension(pt, virtualLength), voropp_container);
    } else {
        throw invalid_argument(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM>
void VoroInterface<DIM>::drawGnuPlot(string fileName) {
    voropp_container->draw_cells_gnuplot(fileName.c_str());
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
inline vector<SingleCell> VoroInterface<DIM>::getSingleCells() {
    if constexpr (DIM == 2) {
        cerr << __PRETTY_FUNCTION__ << endl;
        throw runtime_error("Not programmed yet in dimension 2.");
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
void VoroInterface<DIM>::addInclusion(double* pp, int index,
    voro::voronoicell_neighbor* c,
    vector<smallShape::ConvexPolyhedronInc<DIM>>& polyhedrons) {
    vector<HalfSpace<DIM>> faces = voroInterface_aux::getFaces<DIM>(c); // the faces
    Point<DIM> center; // center of the polyhedron
    Cuboid<DIM> cuboid = voroInterface_aux::getCuboid<DIM>(c, center); // the surrounding cuboid
    if constexpr (DIM == 3) {
        center = { pp[0], pp[1], pp[2] };
    } else  if constexpr (DIM == 2) {
        center = { pp[0], pp[1] };
    }
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
    Point<DIM> center; // center of the polyhedron
    if constexpr (DIM == 3) {
        center = { pp[0], pp[1], pp[2] };
    } else  if constexpr (DIM == 2) {
        center = { pp[0], pp[1] };
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
Cuboid<DIM> voroInterface_aux::getCuboid(voro::voronoicell_neighbor* cell, const Point<DIM>& center) {
    if (not cell) {
        return Cuboid<DIM>(center, center);
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
        int index = voropp_container->id[ijk][q]; // index of the polyhedron
        pp = voropp_container->p[ijk] + voropp_container->ps * q;
        if (voropp_container->compute_cell(c, vl)) {
            my_function(pp, index, &c);
        } else {
            my_function(pp, index, nullptr);
        }
    } while (vl.inc());
}

} // namespace voroInterface
} // namespace merope


#endif /* VOROINTERFACE_IXX_ */
