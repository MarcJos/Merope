//! Copyright : see license.txt
//!
//! \brief


#include "Voronoi/VoroToMeshGraph.hxx"
#include "../../Geometry/include/GeomTools_1.hxx"
#include "Mesh/GeoObjects.hxx"
#include "Mesh/MeshStructure.hxx"


#include "container.hh"
#include "pre_container.hh"


namespace merope {
namespace microToMesh {
namespace voroTranslater {
using namespace mesh::geoObjects;

VoroToMeshGraph::VoroToMeshGraph(array<double, 3> L, const vector<Sphere<DIM>>& centerTessels_, double adimensionnalMergeDistance_) : InsideTorus<3>(L), adimensionnalMergeDistance{ adimensionnalMergeDistance_ }, centerTessels(centerTessels_), rawMeshGraph{} {
    auxi_compute();
}

void VoroToMeshGraph::auxi_compute() {
    // dilate for small/large L
    double dilation = 1. / pow(auxi_function::productOf<double>(getL()), 1. / 3);
    double antidilation = 1. / dilation;
    auto theSpheres = centerTessels;
    for (auto& sphere : theSpheres) {
        sphere.radius *= dilation;
        sphere.center *= dilation;
    }
    auto dilated_L = dilation * getL();
    voroInterface::VoroInterface<3> voroInt(dilated_L, theSpheres);
    auto cells = voroInt.getSingleCells();
    // dilate back
    for (auto& c : cells) {
        c.linTransform(antidilation);
    }
    rawMeshGraph = microToMesh::voroTranslater::auxi::computeFromListOfCells(getL(), cells, adimensionnalMergeDistance);
}

mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const voroInterface::SingleCell& singleCell) {
    mesh::meshStructure::VoroMesh_UnStructureData<3> result{ create_array<3>(0.), {}, {}, {}, {}, {}, {}, {}, {} };
    for (size_t i = 0; i < singleCell.vertices.size(); i++) {
        const auto& pt = singleCell.vertices[i];
        result.vecPoint.push_back(mesh::geoObjects::GeoPoint<3>(i, pt));
    }
    long k = 1;
    vector<Identifier> leaves{};
    vector<Identifier> surfaceLoopLeaves{};
    for (size_t index_surf = 0; index_surf < singleCell.faceVertices.size(); index_surf++) {
        const auto& vertices = singleCell.faceVertices[index_surf];
        vector<Identifier> curveLoopLeaves{};
        for (size_t j = 0; j < vertices.size(); j++) {
            size_t j_ = (j < vertices.size() - 1) ? j + 1 : 0;
            //
            leaves = { vertices[j], vertices[j_] };
            result.vecEdge.push_back(mesh::geoObjects::Edge(k, leaves));
            curveLoopLeaves.push_back(k);
            k++;
            //
        }
        result.vecCurveLoop.push_back(mesh::geoObjects::CurveLoop(index_surf + 1, curveLoopLeaves));
        //
        leaves = { static_cast<Identifier>(index_surf + 1) };
        result.vecSurface.push_back(mesh::geoObjects::Surface(index_surf + 1, leaves));
        surfaceLoopLeaves.push_back(index_surf + 1);
        //
    }
    result.vecSurfaceLoop.push_back(mesh::geoObjects::SurfaceLoop(1, surfaceLoopLeaves));
    //
    leaves = { 1 };
    result.vecSolid.push_back(mesh::geoObjects::Solid(1, leaves));
    //
    return result;
}

mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const Point<3>& center, const vector<HalfSpace<3> >& center_to_faces, const Point<3>& L) {
    voro::voronoicell_neighbor voroCell{};
    voroCell.init(-L[0] * 0.5, L[0] * 0.5, -L[1] * 0.5, L[1] * 0.5, -L[2] * 0.5, L[2] * 0.5);
    // build the cell
    for (const auto& face : center_to_faces) {
        voroCell.plane(face.vec()[0], face.vec()[1], face.vec()[2], 2 * face.c());
    }
    //
    auto singleCell = voroInterface::SingleCell(1, center, &voroCell);
    return getRawMeshGraph(singleCell);
}

mesh::meshStructure::VoroMesh_UnStructureData<3> auxi::computeFromListOfCells(
    Point<3> L, const vector<merope::voroInterface::SingleCell>& outputVoroPP, double adimensionnalMergeDistance) {
    constexpr unsigned short DIM = 3;
    AmbiantSpace::Tore<DIM> torus(L);
    double epsilon_0 = adimensionnalMergeDistance * pow(auxi_function::productOf<double>(L), 1. / DIM);

    std::map<Identifier, size_t> dictSingleCells{};
    //! store the link between identifier a of SingleCell -> index of outputVoroPP
    for (size_t i = 0; i < outputVoroPP.size(); i++) {
        dictSingleCells[outputVoroPP[i].identifier] = i;
    }

    vector<Point<DIM>> vertices{};
    //! store the list of all vertices;
    std::map<Identifier, vector<size_t>> dictPoints{};
    //! store the link between indexes of points and vertices in a SingleCell
    size_t index_point = 0;
    for (const auto& singleCell : outputVoroPP) {
        dictPoints[singleCell.identifier] = {};
        for (const auto& point : singleCell.vertices) {
            vertices.push_back(point);
            dictPoints[singleCell.identifier].push_back(index_point);
            index_point++;
        }
    }

    vector<pair<size_t, size_t>> points_replacementList = {}; // id_pt1, id_pt2
    //! store the replacement list of points
    // compute the replacement list
    vector<pair<pair<Identifier, size_t>, pair<Identifier, size_t>>> correspondingSurfaces{};
    for (const auto& singleCell_0 : outputVoroPP) {
        Identifier id_singleCell_0 = singleCell_0.identifier;
        for (size_t face_id_0 = 0; face_id_0 < singleCell_0.faceVertices.size(); face_id_0++) {
            Identifier id_singleCell_1 = abs(singleCell_0.neighbors[face_id_0]);
            auto& singleCell_1 = outputVoroPP[dictSingleCells[id_singleCell_1]];
            //
            vector<size_t> candidatesFaces = voroInterface::auxi::correspondingFaces<false>([&](Identifier id) {
                if (dictSingleCells.find(abs(id)) != dictSingleCells.end()) {
                    return &(outputVoroPP[dictSingleCells[abs(id)]]);
                } else {
                    return static_cast<const merope::voroInterface::SingleCell*>(nullptr);
                }
                },
                id_singleCell_0, face_id_0).second;

            if (candidatesFaces.size() == 0) {
                std::cerr << "An entire face should be merged" << endl;
                std::cerr << "warning : correct checks are not implemented" << endl;
                //! merge all points of the face
                for (size_t ii = 1; ii < singleCell_0.faceVertices[face_id_0].size(); ii++) {
                    points_replacementList.push_back({ dictPoints[id_singleCell_0][singleCell_0.faceVertices[face_id_0][ii]], dictPoints[id_singleCell_0][singleCell_0.faceVertices[face_id_0][0]] });
                }
            } else if (candidatesFaces.size() == 1) {
                size_t face_id_1 = candidatesFaces[0];
                correspondingSurfaces.push_back({ {id_singleCell_0, face_id_0}, {id_singleCell_1, face_id_1} });
                if (abs(singleCell_1.neighbors[face_id_1]) != singleCell_0.identifier) {
                    throw runtime_error("Face index error. Unexpected!");
                }
                if (singleCell_0.faceVertices[face_id_0].size()
                    != singleCell_1.faceVertices[face_id_1].size()) {
                    cerr << "Different number of points in each face.";
                    Merope_error_not_done();
                }

                // remove too close points in each face
                vector<Identifier> id_points_0 = {}, id_points_1 = {};

                for (size_t iii = 0; iii < singleCell_0.faceVertices[face_id_0].size(); iii++) {
                    id_points_0.push_back(dictPoints[id_singleCell_0][singleCell_0.faceVertices[face_id_0][iii]]);
                    id_points_1.push_back(dictPoints[id_singleCell_1][singleCell_1.faceVertices[face_id_1][iii]]);
                }

                // merge close points
                {
                    auto merge_points = [&](auto& id_points, double epsilon_0_) {
                        for (long i = 0; i < id_points.size(); i++) {
                            size_t j = auxi_function::fast_modulo(i + 1, id_points_0.size()); // i+1 with periodicity
                            size_t id_pt_i = id_points[i];
                            size_t id_pt_j = id_points[j];
                            if (geomTools::normeCarre<DIM>(vertices[id_pt_i] - vertices[id_pt_j]) < epsilon_0_ * epsilon_0_) {
                                points_replacementList.push_back({ id_pt_i , id_pt_j });
                            }
                        }
                        };
                    merge_points(id_points_0, epsilon_0);
                    merge_points(id_points_1, epsilon_0);
                }


                long i_begin = 0;
                {
                    double energy = 0.5 * numeric_limits<double>::max();
                    for (long i = 0; i < id_points_0.size(); i++) {
                        double local_energy = 0;
                        for (long j = 0; j < id_points_1.size(); j++) {
                            size_t k = auxi_function::fast_modulo(i - j, id_points_0.size());
                            const auto& pt_0 = vertices[id_points_0[j]];
                            const auto& pt_1 = vertices[id_points_1[k]];
                            local_energy += torus.distanceCarre(pt_0, pt_1);
                        }
                        if (local_energy < energy) {
                            energy = local_energy;
                            i_begin = i;
                        }
                    }
                }
                for (long j = 0; j < id_points_1.size(); j++) {
                    size_t k = auxi_function::fast_modulo(i_begin - j, id_points_0.size());
                    points_replacementList.push_back({ id_points_0[j], id_points_1[k] });
                }
            } else { // if (candidatesFaces.size() > 1)
                Merope_assert(false, "More than one candidate face!");
            }
        }
    }

    // build merge list
    // point_id_0 -> point_id_1 where point_0 is to be merged into point_id_1 (in the torus)
    vector<unordered_set<size_t>> pointsEquivalenceClasses = merope::cppFunctions::computeEquivalenceClass<size_t>(points_replacementList);
    std::map<size_t, size_t> merge_list = {};
    for (const auto& setPtsId : pointsEquivalenceClasses) {
        auto minimum = *(std::min_element(setPtsId.begin(), setPtsId.end()));
        for (const auto& elem : setPtsId) {
            if (elem != minimum) {
                merge_list[elem] = minimum;
            }
        }
    }


    for (const auto& pair : merge_list) {
        size_t id_replacement = get<1>(pair);
        Point<DIM> replacingPoint = vertices[id_replacement];
        // for periodicity
        auto deltaL = vertices[get<0>(pair)] - replacingPoint;
        for (size_t d = 0; d < DIM; d++) {
            deltaL[d] = std::round(deltaL[d] / L[d]) * L[d];
        }
        replacingPoint += deltaL;
        vertices[get<0>(pair)] = replacingPoint;
    }


    //!
    mesh::meshStructure::VoroMesh_UnStructureData<3> result{};
    result.L = L;

    for (Identifier id = 0; id < vertices.size(); id++) {
        result.vecPoint.push_back(GeoPoint<3>(id + 1, vertices[id]));
    }

    //
    std::map<size_t, vector<size_t>> periodicCopies{};
    std::map<pair<Identifier, size_t>, size_t> surfPair_to_id{};
    for (const auto& singleCell : outputVoroPP) {
        const auto& idPoints = dictPoints[singleCell.identifier];
        vector<Identifier> realPointsId = {};
        for (auto id : idPoints) {
            size_t new_id = id;
            auto found = merge_list.find(id);
            if (found != merge_list.end()) {
                if (geomTools::normeCarre<3>(vertices[found->second] - vertices[id]) < 0.01 * geomTools::normeCarre<3>(L)) {
                    // no periodicity
                    new_id = found->second;
                } else {
                    // periodicity
                    auto found_periodic = periodicCopies.find(found->second);
                    if (found_periodic == periodicCopies.end()) {
                        periodicCopies[found->second] = { new_id };
                    } else {
                        bool success = false;
                        for (auto other_id : found_periodic->second) {
                            if (geomTools::normeCarre<3>(vertices[new_id] - vertices[other_id]) < 0.01 * geomTools::normeCarre<DIM>(L)) {
                                new_id = other_id;
                                success = true;
                                break;
                            }
                        }
                        if (not success) {
                            periodicCopies[found->second].push_back(new_id);
                        }
                    }
                }
            }
            realPointsId.push_back(new_id + 1);
        }
        //
        vector<Identifier> surfaces = {};
        for (size_t id_face_in_cell = 0; id_face_in_cell < singleCell.faceVertices.size(); id_face_in_cell++) {
            const auto& fv = singleCell.faceVertices[id_face_in_cell];
            vector<Identifier> edges = {};
            for (size_t i = 0; i < fv.size(); i++) {
                size_t j = auxi_function::fast_modulo(i + 1, fv.size());
                Identifier id_edge = result.vecEdge.size() + 1;
                edges.push_back(id_edge);
                result.vecEdge.push_back(Edge(id_edge, { realPointsId[fv[i]], realPointsId[fv[j]] }));
            }
            Identifier id_curveLoop = result.vecCurveLoop.size() + 1;
            result.vecCurveLoop.push_back(CurveLoop(id_curveLoop, edges));
            Identifier id_surface = result.vecSurface.size() + 1;
            result.vecSurface.push_back(Surface(id_surface, { id_curveLoop }));
            surfPair_to_id[{ singleCell.identifier, id_face_in_cell }] = id_surface;
            surfaces.push_back(id_surface);
        }
        Identifier id_surfaceLoop = singleCell.identifier + 1;
        result.vecSurfaceLoop.push_back(SurfaceLoop(id_surfaceLoop, surfaces));
        Identifier id_solid = singleCell.identifier + 1;
        result.vecSolid.push_back(Solid(id_solid, { id_surfaceLoop }));
    }

    for (const auto& pair : periodicCopies) {
        Identifier identifier = pair.first + 1;
        vector<Identifier> leaves = { identifier };
        for (const auto& id : pair.second) {
            leaves.push_back(id + 1);
        }
        result.vecPerPoint.push_back(PerPoint(identifier, leaves));
    }


    /*
    // what to do with that info?
    ofstream file("Pairs.txt");
    for (const auto& perSurf : correspondingSurfaces) {
        file << "Persurface " << surfPair_to_id[get<0>(perSurf)] << " : " << surfPair_to_id[get<1>(perSurf)] << endl;

    }
    */
    return result;
}

}  // namespace voroTranslater
}  // namespace microToMesh
}  // namespace merope

