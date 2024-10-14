//! Copyright : see license.txt
//!
//! \brief
//

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "MicroToMesh/MeshGenerator.hxx"
#include "MicroToMesh/SphereMesh.hxx"
#include "Mesh/GmshWriter.hxx"
#include "Voronoi/VoroToMeshGraph.hxx"


#include "MeropeNamespace.hxx"


namespace merope {
namespace mesh {
namespace generator {

void MeshGenerator::write(string nameFile) const {
    ofstream file(nameFile);
    this->write(file);
}

void MeshGenerator::write(std::ostream& f) const {
    auto geoPerStructure = this->meshFrom_MultiInclusions();
    geoPerStructure.isCoherent();
    //
    auxi::create_physical_surfaces(geoPerStructure);
    auxi::remove_ignored_phase(geoPerStructure, ignore_interior);
    auxi::check_phases(geoPerStructure);
    //
    gmsh_writer::auxi::writePreamble(f);
    auxi::write_all(f, geoPerStructure);
    gmsh_writer::auxi::writeEnd(f, this->get_nameOutput(), this->getMeshSize(), this->getMeshOrder(), this->binaryOutput);
}

std::tuple<mesh::meshStructure::VoroMesh_UnStructureData<3>,
    std::map<PhaseType, mesh::geoObjects::PhysicalVolume>,
    vector<Identifier>> MeshGenerator::computeRawDataInclusions() const {
    //
    std::map<PhaseType, geoObjects::PhysicalVolume> physicalVolumes{};
    //
    mesh::meshStructure::VoroMesh_UnStructureData<3> rawMeshData{};
    rawMeshData.L = this->multiInclusions->getL();
    //
    vector<Identifier> outerSurfaces{};
    //// Prepare for the loop
    Identifier current_identifier = 1;
    //
    auto my_func = [&](const auto& polyhedron) {
        long nbInclusions = polyhedron.getInnerInclusions().size();
        Identifier last_non_empty_inner_surfaceLoop = -1;
        for (auto i = nbInclusions - 1; i >= 0; i--) {
            auto polyhedronMeshData = auxi::getRawMeshGraph(polyhedron, rawMeshData.L, i);
            if (polyhedronMeshData.vecSolid.size() > 0) {  // the inclusion is of volume > 0
                polyhedronMeshData.shiftIndices(current_identifier);
                if (last_non_empty_inner_surfaceLoop > 0) {
                    polyhedronMeshData.vecSolid[0].leaves.push_back(last_non_empty_inner_surfaceLoop);
                }
                rawMeshData.append(polyhedronMeshData);
                { // update PhysicalVolumes
                    auto phase = polyhedron.getPhaseForVoxellation(i);
                    auto id = polyhedronMeshData.vecSolid[0].identifier;
                    if (physicalVolumes.find(phase) != physicalVolumes.end()) {
                        physicalVolumes.at(phase).leaves.push_back(id);
                    } else {
                        physicalVolumes.insert(make_pair(phase, geoObjects::PhysicalVolume(phase, { id })));
                    }
                }
                last_non_empty_inner_surfaceLoop = polyhedronMeshData.vecSolid[0].leaves[0];
                // update outer surfaces
                if (i == 0) outerSurfaces.push_back(polyhedronMeshData.vecSurfaceLoop[0].identifier); // outer surface
                // update indices
                current_identifier = polyhedronMeshData.getMaxIndex() + 1;
            }
        }
        };
    //
    this->multiInclusions->apply_on_all([&](const auto& vecObj) {
        for (const auto& obj : vecObj) {
            my_func(obj);
        }
        });
    //
    return make_tuple(rawMeshData, physicalVolumes, outerSurfaces);
}


mesh::meshStructure::VoroMesh_Periodic_Physical<3> MeshGenerator::meshFrom_LaguerreTess() const {
    auto phy_raw = computeRawDataInclusions();
    const auto& rawMeshData = get<0>(phy_raw);
    const auto& physicalVolumes = get<1>(phy_raw);
    //
    mesh::meshStructure::VoroMesh_Periodic_Physical<3>
        geoPerStructure(mesh::meshStructure::VoroMesh_Periodic<3>(rawMeshData, adimensionnalMergeDistance_0, adimensionnalMergeDistance_1));
    geoPerStructure.dictPhysicalVolume = physicalVolumes;
    return geoPerStructure;
}

mesh::meshStructure::VoroMesh_Periodic_Physical<3> MeshGenerator::getBoundaryMesh() const {
    vector<Sphere<3>> theSpheres{};
    //
    if (this->multiInclusions->template get<smallShape::SphereInc<3>>().size() > 0) {
        const auto& sphereInc = this->multiInclusions->template get<smallShape::SphereInc<3>>();
        theSpheres.resize(sphereInc.size());
        for (size_t i = 0, l = sphereInc.size(); i < l; i++) {
            theSpheres[i] = sphereInc[i].getInnerInclusions()[0];
        }
    } else if (this->multiInclusions->template get<smallShape::CylinderInc<3>>().size() > 0) {
        { // MESSAGE
            for (size_t i = 0; i < 10; i++) {
                cerr << "WARNING" << endl;
            }
            cerr << "Implicitly assume that the cylinders belong to non-intersecting spheres!" << endl;
            for (size_t i = 0; i < 10; i++) {
                cerr << "WARNING" << endl;
            }
        }
        for (const auto& cylinderInc : this->multiInclusions->template get<smallShape::CylinderInc<3>>()) {
            const auto& outer_cylinder = cylinderInc.getInnerInclusions()[0];
            double half_diagonal = sqrt(0.25 * geomTools::normeCarre<3>(outer_cylinder.axis[1] - outer_cylinder.axis[0])
                + outer_cylinder.radius * outer_cylinder.radius);
            theSpheres.push_back(Sphere<3>(outer_cylinder.center(), half_diagonal, 0));
        }
    } else {
        throw runtime_error("No known boundary for such a shape");
    }
    //
    microToMesh::voroTranslater::VoroToMeshGraph voroTranslater(multiInclusions->getL(), theSpheres);
    mesh::meshStructure::VoroMesh_UnStructureData<3> rawMeshData = voroTranslater.getMeshData();
    mesh::meshStructure::VoroMesh_Periodic_Physical<3>
        geoPerStructure(mesh::meshStructure::VoroMesh_Periodic<3>(rawMeshData, adimensionnalMergeDistance_0, adimensionnalMergeDistance_1));
    geoPerStructure.restrictEnveloppe();
    geoPerStructure.isStronglyCoherent();
    return geoPerStructure;
}

mesh::meshStructure::VoroMesh_Periodic_Physical<3> MeshGenerator::meshFrom_MultiInclusions() const {
    if (multiInclusions->template get<smallShape::ConvexPolyhedronInc<3>>().size() > 0) {
        if ((multiInclusions->is_there_matrix())) {
            cerr << "WARNING" << endl; // warning
            cerr << __PRETTY_FUNCTION__ << endl;
            cerr << "I expect either a laguerre tessellation without a matrix or a matrix with inclusions" << endl;
        }
        return this->meshFrom_LaguerreTess();
    } else {
        return this->meshFrom_MatrixInclusions();
    }
}

mesh::meshStructure::VoroMesh_Periodic_Physical<3> MeshGenerator::meshFrom_MatrixInclusions() const {
    mesh::meshStructure::VoroMesh_Periodic_Physical<3> geoPerStructure = getBoundaryMesh();
    //!
    auto phy_raw = computeRawDataInclusions();
    auto& rawMeshData = get<0>(phy_raw);
    auto& physicalVolumes = get<1>(phy_raw);
    auto& inclusionSurfacesLoop = get<2>(phy_raw);
    //
    //!
    Identifier Last_Identifier = geoPerStructure.getMaxIndex() + 1;
    { // shift everyone by Last_Identifier
        rawMeshData.shiftIndices(Last_Identifier);
        for (auto& [id, phyV] : physicalVolumes) {
            for (auto& id_solid : phyV.leaves) {
                id_solid += Last_Identifier;
            }
        }
        for (auto& id : inclusionSurfacesLoop) {
            id += Last_Identifier;
        }
    }

    // matrix
    if (this->multiInclusions->is_there_matrix()) {
        Identifier matrixPhase = this->multiInclusions->getMatrixPhase();
        physicalVolumes.insert(make_pair(matrixPhase, geoObjects::PhysicalVolume(matrixPhase, { 1 })));
        inclusionSurfacesLoop.insert(inclusionSurfacesLoop.begin(), geoPerStructure.dictSolid.at(1).leaves[0]);
        geoPerStructure.dictSolid.at(1).leaves = inclusionSurfacesLoop;
    }
    //!
    geoPerStructure.dictPhysicalVolume = physicalVolumes;
    geoPerStructure.add_raw_mesh_data(rawMeshData);
    return geoPerStructure;
}

namespace auxi {

void write_all(std::ostream& f, const mesh::meshStructure::VoroMesh_Periodic_Physical<3>& geoPerStructure) {
    gmsh_writer::auxi::writeDict("Points", geoPerStructure.dictPoint, f);
    gmsh_writer::auxi::writeDict("Edges", geoPerStructure.dictEdge, f);
    gmsh_writer::auxi::writeDict("CurveLoops", geoPerStructure.dictCurveLoop, f);
    gmsh_writer::auxi::writeDict("Surfaces", geoPerStructure.dictSurface, f);
    gmsh_writer::auxi::writeDict("SurfaceLoop", geoPerStructure.dictSurfaceLoop, f);
    gmsh_writer::auxi::writeDict("Matrix volume", geoPerStructure.dictSolid, f);
    gmsh_writer::auxi::writeDict("Periodic surfaces of the enveloppe", geoPerStructure.dictPerSurface, f);
    gmsh_writer::auxi::writeDict("Physical volumes", geoPerStructure.dictPhysicalVolume, f);
    gmsh_writer::auxi::writeDict("Physical surfaces", geoPerStructure.dictPhysicalSurface, f);
}

void create_physical_surfaces(mesh::meshStructure::VoroMesh_Periodic_Physical<3>& geoPerStructure) {
    Identifier max_phase = 0;
    for (const auto& key_phyVol : geoPerStructure.dictPhysicalVolume) {
        auto phase = get<0>(key_phyVol);
        const geoObjects::PhysicalVolume& phyVol = get<1>(key_phyVol);
        geoObjects::PhysicalSurface phySurf(phase, {});
        for (auto id_vol : phyVol.leaves) {
            const auto& solid = geoPerStructure.dictSolid.at(id_vol);
            for (auto id_surfaceLoop : solid.leaves) {
                const auto& surfaceLoop = geoPerStructure.dictSurfaceLoop.at(id_surfaceLoop);
                for (auto id_surf : surfaceLoop.leaves) {
                    phySurf.leaves.push_back(abs(id_surf));
                }
            }
        }
        geoPerStructure.dictPhysicalSurface.insert(make_pair(phase, phySurf));
        max_phase = max(phase, max_phase);
    }
    geoPerStructure.dictPhysicalSurface.insert(make_pair(max_phase + 1, geoPerStructure.getOuterSurface(max_phase + 1)));
}

void check_phases(const mesh::meshStructure::VoroMesh_Periodic_Physical<3>& geoPerStructure) {
    //! check no phase at 0
    for (const auto& key_fv : geoPerStructure.dictPhysicalVolume) {
        if (get<0>(key_fv) < 1) {
            cerr << __PRETTY_FUNCTION__ << endl;
            cerr << "I see a phase " << get<0>(key_fv) << endl;
            cerr << "The phases are forced to be strictly larger than 0" << endl;
            throw runtime_error("Incorrect phase!");
        }
    }
}

void remove_ignored_phase(mesh::meshStructure::VoroMesh_Periodic_Physical<3>& geoPerStructure, const std::map<PhaseType, bool>& ignore_interior) {
    //! remove from dict_solid
    for (const auto key_fv : geoPerStructure.dictPhysicalVolume) {
        if (ignore_interior.find(key_fv.first) != ignore_interior.end()) {
            for (auto s : key_fv.second.leaves) {
                geoPerStructure.dictSolid.erase(s);
            }
        }
    }
    //! remove from physicalVolumes
    std::erase_if(geoPerStructure.dictPhysicalVolume, [&ignore_interior](const auto& c) {
        return ignore_interior.find(c.first) != ignore_interior.end();
        });
}

} // namespace  auxi

}  // namespace generator
}  // namespace mesh
}  // namespace merope
