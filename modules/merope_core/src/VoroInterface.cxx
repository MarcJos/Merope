//! Copyright : see license.txt
//!
//! \brief

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "Geometry/GeomTools.hxx"

#if defined(_WIN32) || defined(WIN32) // ugly for Eclipse
#include "../../voro-plus-plus/src/container.hh"
#include "../../voro-plus-plus/src/pre_container.hh"
#else
#include "container.hh"
#include "pre_container.hh"
#endif

#include "Voronoi/VoroInterface.hxx"


#include "MeropeNamespace.hxx"


namespace merope {
namespace voroInterface {

void voroInterface_aux::buildTessellation(
    const vector<Sphere<3>>& centerTessels, array<double, 3> L,
    voro::container_poly*& voropp_container, array<bool, 3> periodicity) {
    voro::pre_container_poly precont(0., L[0], 0., L[1], 0., L[2],
        periodicity[0], periodicity[1], periodicity[2]);
    // Place the centers
    for (size_t i = 0; i < centerTessels.size(); i++) {
        const auto& sphere = centerTessels[i];
        precont.put(i, sphere.center[0], sphere.center[1], sphere.center[2],
            sphere.radius);
    }
    // find the nx, ny, nz
    int n_x, n_y, n_z;
    precont.guess_optimal(n_x, n_y, n_z);
    // Computes the tessellation
    voropp_container = new voro::container_poly(0., L[0], 0., L[1], 0., L[2],
        n_x, n_y, n_z, periodicity[0], periodicity[1], periodicity[2],
        VOROPP_INIT_MEM);
    precont.setup(*voropp_container);
    voropp_container->compute_all_cells();
}

vector<Point<3> > voroInterface_aux::breakList(const vector<double>& v) {
    size_t s = v.size();
    size_t imax = s / 3;
    if (s != 3 * imax) {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
    vector<Point<3>> result(imax);
    for (size_t i = 0; i < imax; i++) {
        result[i] = Point<3>{ v[3 * i], v[3 * i + 1], v[3 * i + 2] };
    }
    return result;
}

vector<HalfSpace<3>> voroInterface_aux::getFaces3D(voro::voronoicell_neighbor* cell) {
    SingleCell singleCell(0, create_array<3>(0.), cell);
    auto numFaces = singleCell.faceNormal.size();
    // prepare the faces
    vector<HalfSpace<3>> faces{};
    faces.reserve(numFaces);
    //////////////////////////////////////////////////////////////////////////////
    for (size_t j = 0; j < numFaces; j++) {
        {
            const Point<3>& pointOnFace = singleCell.vertices.at(singleCell.faceVertices[j][0]);
            const Point<3>& normal = singleCell.faceNormal[j];
            faces.push_back(HalfSpace<3>(normal, pointOnFace));
        }
    }
    return faces;
}

vector<HalfSpace<2>> voroInterface_aux::getFaces2D(
    voro::voronoicell_neighbor* cell) {
    vector<HalfSpace<3>> halfSpaces = getFaces3D(cell);
    vector<HalfSpace<2>> result{};
    for (const auto& hf : halfSpaces) {
        Point<2> newNormal = { hf.vec()[0], hf.vec()[1] };
        if (geomTools::normeCarre<2>(newNormal) > 0.9 * geomTools::normeCarre<3>(hf.vec())) {
            result.push_back(HalfSpace<2>(newNormal, hf.c()));
        }
    }
    return result;
}

vector<Point<3> > voroInterface_aux::getRenormalizedVertices(
    voro::voronoicell_neighbor* cell) {
    vector<double> vecTransfer{ };
    vector<int> neigh;
    cell->neighbors(neigh);
    cell->vertices(vecTransfer);
    return voroInterface_aux::breakList(vecTransfer);
}

vector<Point<3>> voroInterface_aux::getNormals(voro::voronoicell_neighbor* cell) {
    vector<double> vecTransfer{ }; //used to get the information from voro++
    cell->normals(vecTransfer);
    return voroInterface_aux::breakList(vecTransfer);
}

Point<3> voroInterface_aux::extendDimension(const Point<2>& oldPoint, const array<double, 3>& L) {
    return Point<3>{oldPoint[0], oldPoint[1], 0.5 * L[2]};
}

vector<Point<3> > voroInterface_aux::getVertices(voro::voronoicell_neighbor* cell, Point<3> center) {
    auto result = voroInterface_aux::getRenormalizedVertices(cell);
    for (auto& pt : result) {
        pt = pt + center;
    }
    return result;
}

int voroInterface_aux::findTessel(const Point<3>& pt,
    voro::container_poly* voropp_container) {
    double a, b, c;
    int indexTessel;
    bool success = voropp_container->find_voronoi_cell(pt[0], pt[1], pt[2], a,
        b, c, indexTessel);
    if (success) {
        return indexTessel;
    } else {
        throw runtime_error(__PRETTY_FUNCTION__);
    }
}

vector<vector<Identifier>> voroInterface_aux::getFacesToVertices(
    voro::voronoicell_neighbor* cell) {
    //! extract the information
    vector<int> face2vertices{};
    cell->face_vertices(face2vertices);
    //! convert the information
    vector<vector<Identifier>> result{};
    size_t i = 0; // index of the face
    size_t k = 0; // position of the reading head in face2vertices
    while (k < face2vertices.size()) {
        size_t jmax = face2vertices[k]; k++; // read the nb of vertices in the face
        result.push_back(vector<Identifier>(jmax));
        for (size_t j = 0; j < jmax; j++) {
            result[i][j] = face2vertices[k]; k++;
        }
        i++;
    }
    return result;
}

SingleCell::SingleCell(Identifier identifier_, const Point<3>& center_,
    voro::voronoicell_neighbor* c) :
    SingleCell{ identifier_, center_,
            voroInterface_aux::getVertices(c, center_),
            voroInterface_aux::getNormals(c),
            voroInterface_aux::getFacesToVertices(c),
            voroInterface_aux::getNeighbors(c) } {}

// SingleCell
void SingleCell::print(ostream& os) const {
    os << "Identifier : " << this->identifier << endl;
    os << "Center : "; auxi_function::writeVectorToString(this->center, os); os << endl;
    os << "Vertices : ";
    for (const auto& v : this->vertices) {
        auxi_function::writeVectorToString(v, os); os << " ; ";
    } os << endl;
    os << "faceNormal : ";
    for (const auto& v : this->faceNormal) {
        auxi_function::writeVectorToString(v, os); os << " ; ";
    } os << endl;
    os << "faceVertices : ";
    for (const auto& v : this->faceVertices) {
        auxi_function::writeVectorToString(v, os); os << " ; ";
    } os << endl;
    os << "neighbords : "; auxi_function::writeVectorToString(this->neighbors, os); os << endl;
}

void SingleCell::correctNormals() {
    for (size_t i = 0; i < faceNormal.size(); i++) {
        if (geomTools::normeCarre<3>(faceNormal[i]) < 1e-10) {
            // fixme
            /*
            faceNormal[i] = geomTools::prodVec<3>(vertices[faceVertices[i][0]] - vertices[faceVertices[i][1]], vertices[faceVertices[i][1]] - vertices[faceVertices[i][2]]) ;
             */
            faceNormal[i] = vertices[faceVertices[i][0]] - center;
            geomTools::renormalize<3>(faceNormal[i]);
        }
    }
}

vector<Identifier> voroInterface_aux::getNeighbors(voro::voronoicell_neighbor* cell) {
    vector<int> temp_result{};
    cell->neighbors(temp_result);
    vector<Identifier> result{};
    std::copy(temp_result.begin(), temp_result.end(), std::back_inserter(result));
    return result;
}

SingleCell::SingleCell(Identifier identifier_, const Point<3>& center_,
    const vector<Point<3> >& vertices_,
    const vector<Point<3> >& faceNormal_,
    const vector<vector<Identifier> >& faceVertices_,
    const vector<Identifier>& neighbors_) :
    identifier{ identifier_ }, center(center_), vertices(vertices_), faceNormal(faceNormal_), faceVertices(faceVertices_),
    neighbors(neighbors_) {
    correctNormals(); // for avoiding 0 normals
}

} // namespace voroInterface
} // namespace merope


