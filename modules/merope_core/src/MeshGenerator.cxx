//! Copyright : see license.txt
//!
//! \brief
//

#include "../../AlgoPacking/src/StdHeaders.hxx"

#include "MicroToMesh/MeshGenerator.hxx"
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
    if (multiInclusions->getPolyhedrons().size() == 0) {
        this->writeSphericalInclusions(f);
    }
    else {
        this->writeLaguerreTess(f);
    }
    gmsh_writer::auxi::writeEnd(f, this->getMeshSize(), this->getMeshOrder(), this->binaryOutput);
}

void MeshGenerator::writeLaguerreTess(std::ostream& f) const {
    std::map<PhaseType, vector<Identifier>> physicalVolumes_id{};
    //!
    mesh::meshStructure::VoroMesh_UnStructureData<3> rawMeshData{};
    rawMeshData.L = this->multiInclusions->getL();
    //// Prepare for the loop
    Identifier current_identifier = 1;
    //
    for (const auto& polyhedron : this->multiInclusions->getPolyhedrons()) {
        long nbInclusions = polyhedron.getInnerInclusions().size();
        Identifier last_non_empty_inner_surfaceLoop = -1;
        for (auto i = nbInclusions - 1; i >= 0; i--) {
            auto polyhedronMeshData = microToMesh::voroTranslater::getRawMeshGraph(polyhedron.center, polyhedron.getInnerInclusions()[i].faces, rawMeshData.L);
            if (polyhedronMeshData.vecSolid.size() > 0) { // the inclusion is of volume > 0
                polyhedronMeshData.shiftIndices(current_identifier);
                if (last_non_empty_inner_surfaceLoop > 0) {
                    polyhedronMeshData.vecSolid[0].leaves.push_back(last_non_empty_inner_surfaceLoop);
                }
                rawMeshData.append(polyhedronMeshData);
                physicalVolumes_id[polyhedron.getPhaseForVoxellation(i)].push_back(polyhedronMeshData.vecSolid[0].identifier);
                // update indices
                last_non_empty_inner_surfaceLoop = polyhedronMeshData.vecSolid[0].leaves[0];
                current_identifier = polyhedronMeshData.getMaxIndex() + 1;
            }
        }
    }
    //
    mesh::meshStructure::VoroMesh_Periodic<3> voroMeshPer(rawMeshData, adimensionnalMergeDistance_0, adimensionnalMergeDistance_1);
    std::map<PhaseType, geoObjects::PhysicalVolume> physicalVolumes{};
    for (const auto& [phase, ids] : physicalVolumes_id) {
        bool nonVoidSurface = false;
        for (auto id : ids) {
            // fixme : unefficient sort
            if (voroMeshPer.dictSolid.find(id) != voroMeshPer.dictSolid.end()) {
                if (nonVoidSurface) {
                    physicalVolumes.at(phase).leaves.push_back(id);
                }
                else {
                    physicalVolumes.insert(make_pair(phase, geoObjects::PhysicalVolume(phase, { id })));
                    nonVoidSurface = true;
                }
            }
        }
    }
    //
    f.precision(17);
    gmsh_writer::auxi::writePreamble(f);
    //
    gmsh_writer::auxi::writeDict("Points of the enveloppe", voroMeshPer.dictPoint, f);
    gmsh_writer::auxi::writeDict("Edges of the enveloppe", voroMeshPer.dictEdge, f);
    gmsh_writer::auxi::writeDict("CurveLoops of the enveloppe", voroMeshPer.dictCurveLoop, f);
    gmsh_writer::auxi::writeDict("Surfaces of the enveloppe", voroMeshPer.dictSurface, f);
    gmsh_writer::auxi::writeDict("SurfaceLoop of the enveloppe", voroMeshPer.dictSurfaceLoop, f);
    gmsh_writer::auxi::writeDict("Matrix volume", voroMeshPer.dictSolid, f);
    gmsh_writer::auxi::writeDict("Periodic surfaces of the enveloppe", voroMeshPer.dictPerSurface, f);
    gmsh_writer::write(voroMeshPer.getOuterSurface(1), f); // physical outer surface
    //
    gmsh_writer::auxi::writeDict("Physical volumes", physicalVolumes, f);
}

void MeshGenerator::writeSphericalInclusions(std::ostream& f) const {
    // for getting the outer enveloppe
    const auto& sphereInc = this->multiInclusions->getSphereInc();
    vector<Sphere<3>> theSpheres(sphereInc.size());
    for (size_t i = 0, l = sphereInc.size(); i < l; i++) {
        theSpheres[i] = sphereInc[i].getInnerInclusions()[0];
    }
    microToMesh::voroTranslater::VoroToMeshGraph voroTranslater(multiInclusions->getL(), theSpheres);
    mesh::meshStructure::VoroMesh_UnStructureData<3> rawMeshData = voroTranslater.getMeshData();
    mesh::meshStructure::VoroMesh_Periodic<3> geoPerStructure(rawMeshData, adimensionnalMergeDistance_0, adimensionnalMergeDistance_1);
    geoPerStructure.restrictEnveloppe();
    geoPerStructure.isStronglyCoherent();
    //!
    Identifier Last_Identifier = geoPerStructure.getMaxIndex() + 1;
    // get the sphere mesh identifiers for each inclusion
    vector<vector<Sphere<3>>> theSpheres_all(sphereInc.size());
    vector<Identifier> outerSpheres(sphereInc.size());
    std::map<PhaseType, vector<Identifier>> phaseTypeToIdentifier{};
    for (size_t i = 0, l = sphereInc.size(); i < l; i++) {
        theSpheres_all[i] = {};
        for (auto sphere : sphereInc[i].getInnerInclusions()) {
            phaseTypeToIdentifier[sphere.phase].push_back(Last_Identifier);
            sphere.phase = Last_Identifier;
            theSpheres_all[i].push_back(sphere);
            geoPerStructure.dictSolid.at(1).leaves.push_back(Last_Identifier);
            Last_Identifier += mesh::gmsh_writer::NUMBER_MESHCOMPONENT_PER_SPHERE;
        }
    }
    // write everyone
    // write the matrix
    f.precision(17);
    gmsh_writer::auxi::writePreamble(f);
    // write the spheres
    for (const auto& listSpheres : theSpheres_all) {
        for (long i = listSpheres.size() - 1; i >= 0; i--) {
            if (i == listSpheres.size() - 1) { // no inner sphere
                gmsh_writer::auxi::write_sphere(listSpheres[i], f);
            }
            else {
                gmsh_writer::auxi::write_sphere(listSpheres[i], f, { listSpheres[i + 1].phase });
            }
        }
    }
    //
    gmsh_writer::auxi::writeDict("Points of the enveloppe", geoPerStructure.dictPoint, f);
    gmsh_writer::auxi::writeDict("Edges of the enveloppe", geoPerStructure.dictEdge, f);
    gmsh_writer::auxi::writeDict("CurveLoops of the enveloppe", geoPerStructure.dictCurveLoop, f);
    gmsh_writer::auxi::writeDict("Surfaces of the enveloppe", geoPerStructure.dictSurface, f);
    gmsh_writer::auxi::writeDict("SurfaceLoop of the enveloppe", geoPerStructure.dictSurfaceLoop, f);
    gmsh_writer::auxi::writeDict("Matrix volume", geoPerStructure.dictSolid, f);
    gmsh_writer::auxi::writeDict("Periodic surfaces of the enveloppe", geoPerStructure.dictPerSurface, f);
    gmsh_writer::write(geoPerStructure.getOuterSurface(1), f); // physical outer surface
    // get physical volumes
    geoObjects::PhysicalVolume  matrix_physicalVolume(1, { 1 });
    gmsh_writer::write(matrix_physicalVolume, f);
    for (const auto& [phase, vect_id] : phaseTypeToIdentifier) {
        gmsh_writer::write(geoObjects::PhysicalVolume(phase, vect_id), f);
    }
}

} // namespace mesh
} // namespace generator
} // namespace merope
