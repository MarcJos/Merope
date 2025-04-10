//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../../../Geometry/include/AmbiantSpace.hxx"

#include "../Mesh/MeshStructure.hxx"
#include "../MultiInclusions/MultiInclusions.hxx"
#include "../Mesh/GmshWriter.hxx"


namespace merope {
namespace mesh {
namespace generator {

class MeshGenerator {
public:
    //! constructor
    MeshGenerator() : meshSize{ 0.05 }, meshOrder{ 1 },
        adimensionnalMergeDistance{ 1.e-5 },
        binaryOutput{ false },
        nameOutput({}),
        ignore_interior{},
        multiInclusions{ nullptr } {}
    //! write the gmsh input into a stream
    //! \param f : outputstream
    void write(std::ostream& f, gmsh_writer::MeshMethod meshMethod = gmsh_writer::MeshMethod::native_gmsh) const;
    //! write the gmsh input into a file
    //! \param nameFile : the name of the file
    void write(string nameFile, gmsh_writer::MeshMethod meshMethod = gmsh_writer::MeshMethod::native_gmsh) const;

    //! @brief : const getter
    const vector<string> get_nameOutput() const { return nameOutput; }
    //! @brief : setter
    void set_nameOutput(vector<string> nameOutput_) { nameOutput = nameOutput_; }
    //! getter
    size_t getMeshOrder() const { return meshOrder; }
    //! setter
    void setMeshOrder(size_t meshOrder_) { this->meshOrder = meshOrder_; }
    //! setter
    void setMultiInclusions(const MultiInclusions<3>& multiInclusions_) { this->multiInclusions = &multiInclusions_; }
    //! getter
    double getMeshSize() const { return meshSize; }
    //! setter
    void setMeshSize(double meshSize_) { this->meshSize = meshSize_; }
    //! setter
    void setAdimMergeDistance(double adimensionnalMergeDistance_) { adimensionnalMergeDistance = adimensionnalMergeDistance_; }
    //! setter
    void setBinaryOutput(bool binaryOutput_) { binaryOutput = binaryOutput_; }
    //! parametrization
    void do_not_mesh(const vector<PhaseType>& phases) { for (const auto id : phases) ignore_interior[id] = true; }

private:
    //! Options
    //! size (=step) of the mesh elements
    double meshSize;
    //! order of the mesh elements
    size_t meshOrder;
    //! triggers the merge distance
    double adimensionnalMergeDistance;
    //! decides whether the output format is binary
    bool binaryOutput;
    //! @brief name of the output
    vector<string> nameOutput;
    //! @brief should I ignore the interior of the surface?
    std::map<PhaseType, bool> ignore_interior;

    //! inner copy of a multiInclusions
    const MultiInclusions<3>* multiInclusions;
private:
    //! @return a mesh structure
    mesh::meshStructure::VoroMesh_Periodic_Physical<3> meshFrom_MultiInclusions() const;
    //! @return a mesh structure
    mesh::meshStructure::VoroMesh_Periodic_Physical<3> meshFrom_MatrixInclusions() const;
    //! @return a mesh structure
    mesh::meshStructure::VoroMesh_Periodic_Physical<3>  meshFrom_LaguerreTess() const;
    //! @return  return a boundary mesh (if necessary)
    mesh::meshStructure::VoroMesh_Periodic_Physical<3> getBoundaryMesh() const;
    //! \return : the raw data of the inclusions
    //! -> data structure containing the inclusions
    //! -> PhysicalVolume of all phases
    //! -> list of outer enveloppes of (outer) inclusions
    std::tuple<mesh::meshStructure::VoroMesh_UnStructureData<3>,
        std::map<PhaseType, mesh::geoObjects::PhysicalVolume>,
        vector<Identifier>> computeRawDataInclusions() const;
};

namespace auxi {
//! @brief checks if no index in geoPerStructure.dictPhysicalVolume is 0
void check_phases(const mesh::meshStructure::VoroMesh_Periodic_Physical<3>& geoPerStructure);
//! @brief remove all the solids corresponding to the indices of physical volumes contained in ignore_interior
//! @param geoPerStructure : mesh informations
//! @param ignore_interior : indices of physicalVolumes to be ignored
void remove_ignored_phase(mesh::meshStructure::VoroMesh_Periodic_Physical<3>& geoPerStructure, const std::map<PhaseType, bool>& ignore_interior);
//! @brief write all components of the mesh
//! @param f : output (usually, text files)
//! @param geoPerStructure : mesh informations
template<gmsh_writer::MeshMethod meshMethod>
void write_all(std::ostream& f, const mesh::meshStructure::VoroMesh_Periodic_Physical<3>& geoPerStructure);
//! @brief from geoPerStructure.dictPhysical volumes, create their boundaries, indicated as physical surfaces.
//! These share the same indices as those.
//! @param geoPerStructure : mesh informations
void create_physical_surfaces(mesh::meshStructure::VoroMesh_Periodic_Physical<3>& geoPerStructure);
//! @brief from a smallShape::SHAPE, for a given layer surface, build the mesh informations
//! @tparam INCLUSION : MicroInclusion_<DIM, OBJ>
//! @param inclusion : inclusion of type smallShape
//! @param L : size of the domain (important for crystals)
//! @param i : index of the layer inside the MicroInclusion_<DIM, OBJ>
template<class INCLUSION>
mesh::meshStructure::VoroMesh_UnStructureData<3> getRawMeshGraph(const INCLUSION& inclusion, const Point<3>& L, size_t i);
} // namespace  auxi

}  // namespace generator
}  // namespace mesh
}  // namespace merope

#include "MeshGenerator.ixx"
