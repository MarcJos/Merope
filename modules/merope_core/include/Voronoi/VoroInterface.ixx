//! Copyright : see license.txt
//!
//! \brief

#ifndef VOROINTERFACE_IXX_
#define VOROINTERFACE_IXX_


#include "../MeropeNamespace.hxx"


namespace merope {
namespace voroInterface {

template<unsigned short DIM>
VoroInterface<DIM>::VoroInterface(array<double, DIM> L,
    const vector<Sphere<DIM>>& centerTessels):
    InsideTorus<DIM>(L), voropp_container{ nullptr }, virtualLength{} {
    setVirtualLength();
    buildTessellation(centerTessels);
}

template<unsigned short DIM>
VoroInterface<DIM>::~VoroInterface() {
    delete voropp_container;
}

template<unsigned short DIM>
int VoroInterface<DIM>::findTessel(const Point<DIM>& pt) {
    if constexpr (DIM == 3) {
        return voroInterface_aux::findTessel(pt, voropp_container);
    }
    else if constexpr (DIM == 2) {
        return voroInterface_aux::findTessel(voroInterface_aux::extendDimension(pt, virtualLength), voropp_container);
    }
    else {
        throw invalid_argument(__PRETTY_FUNCTION__);
    }
}

template<unsigned short DIM>
void VoroInterface<DIM>::drawGnuPlot(string fileName) {
    voropp_container->draw_cells_gnuplot(fileName.c_str());
}

template<unsigned short DIM>
inline void VoroInterface<DIM>::buildTessellation(
    const vector<Sphere<DIM> >& centerTessels) {
    if constexpr (DIM == 3) {
        voroInterface_aux::buildTessellation(centerTessels, this->getL(), voropp_container);
    }
    else if constexpr (DIM == 2) {
        vector<Sphere<3>> newCenterTessels{};
        for (const auto& sph : centerTessels) {
            Point<3> newCenter = voroInterface_aux::extendDimension(sph.center, virtualLength);
            newCenterTessels.push_back(Sphere<3>(newCenter, sph.radius, sph.phase));
        }
        voroInterface_aux::buildTessellation(newCenterTessels, virtualLength, voropp_container, array<bool, 3>{true, true, false}); // no periodicity in the virtual direction (simpler)
    }
}

template<unsigned short DIM>
inline void VoroInterface<DIM>::setVirtualLength() {
    if constexpr (DIM == 2) {
        virtualLength[0] = this->getL()[0];
        virtualLength[1] = this->getL()[1];
        virtualLength[2] = VIRTUAL_LENGTH_MULTIPLIER * min(virtualLength[0], virtualLength[1]);
    }
    else if constexpr (DIM == 3) {
        virtualLength = this->getL();
    }
}

template<unsigned short DIM>
vector<smallShape::ConvexPolyhedronInc<DIM>> VoroInterface<DIM>::getMicroInclusions() {
    vector<smallShape::ConvexPolyhedronInc<DIM>> polyhedrons{ };
    polyhedrons.reserve(voropp_container->total_particles());
    voro::c_loop_all vl(*voropp_container);
    // Heavily inspired from voro::print_custom
    int ijk, q;
    double* pp;
    voro::voronoicell_neighbor c;
    if (vl.start()) do
        if (voropp_container->compute_cell(c, vl)) {
            ijk = vl.ijk;
            q = vl.q;
            int index = voropp_container->id[ijk][q]; // index of the polyhedron
            pp = voropp_container->p[ijk] + voropp_container->ps * q;
            addInclusion(pp, index, &c, polyhedrons);
        } while (vl.inc());
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
    voro::c_loop_all vl(*voropp_container);
    // Heavily inspired from voro::print_custom
    int ijk, q;
    double* pp;
    voro::voronoicell_neighbor c;
    if (vl.start()) do
        if (voropp_container->compute_cell(c, vl)) {
            ijk = vl.ijk;
            q = vl.q;
            int index = voropp_container->id[ijk][q]; // index of the polyhedron
            pp = voropp_container->p[ijk] + voropp_container->ps * q;
            addCell(pp, index, &c, polyhedrons);
        } while (vl.inc());
        return polyhedrons;
}

template<unsigned short DIM>
void VoroInterface<DIM>::addInclusion(double* pp, int index,
    voro::voronoicell_neighbor* c,
    vector<smallShape::ConvexPolyhedronInc<DIM>>& polyhedrons) {
    vector<HalfSpace<DIM>> faces = voroInterface_aux::getFaces<DIM>(c); // the faces
    Cuboid<DIM> cuboid = voroInterface_aux::getCuboid<DIM>(c); // the surrounding cuboid
    Point<DIM> center; // center of the polyhedron
    if constexpr (DIM == 3) {
        center = { pp[0], pp[1], pp[2] };
    }
    else  if constexpr (DIM == 2) {
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
    }
    else  if constexpr (DIM == 2) {
        center = { pp[0], pp[1] };
    }
    polyhedrons.emplace_back(SingleCell(index, center, c));
}



// voroInterface_aux


template<unsigned short DIM>
Cuboid<DIM> voroInterface_aux::getCuboid(voro::voronoicell_neighbor* cell) {
    static_assert(DIM == 2 or DIM == 3);
    Cuboid<3> cuboid3D = smallShape::createCuboid < 3 >(getRenormalizedVertices(cell));
    if constexpr (DIM == 3) {
        return cuboid3D;
    }
    else if constexpr (DIM == 2) {
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

} // namespace voroInterface
} // namespace merope


#endif /* VOROINTERFACE_IXX_ */
