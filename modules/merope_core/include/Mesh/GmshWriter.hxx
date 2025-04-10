//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../Mesh/GeoObjects.hxx"
#include "../Mesh/MeshStructure.hxx"


namespace merope {

namespace mesh {
namespace gmsh_writer {

inline string mesh_size() { return "meshSize"; }
constexpr double MESH_SIZE = 0.05;
constexpr size_t NUMBER_MESHCOMPONENT_PER_SPHERE = 15;

//! Mesh method : either use the native gmsh .geo file, of the package OpenCascade.
//! For the latter, the numbering convention becomes a problem. 
//! Thus, one has to set indices of objects with special commands, such as 
//! p1 = newp; Point(p1) = ...
//! instead of
//! Point(1) = ...
//! Here "newp" selects the correct index for p1.
enum class MeshMethod {
    native_gmsh, OpenCascade
};

//! @return the gmsh name of the object
template<class OBJ>
string name_of(const OBJ& object);
//! @return the gmsh command for new index of the object
template<class OBJ>
string name_new_(const OBJ& object);
//! @return the name letter for the index of the object
template<class OBJ>
string letter_for(const OBJ& object);
//! @return the name letter for the index of the leaves of the object
template<class OBJ>
string letter_for_leaf(const OBJ& object);


//! write a geometrical entity in a gmsh format (.geo)
//! \param object : geometrical entity
//! \param f : output flux
template<MeshMethod meshMethod, class OBJ>
void write(const OBJ& object, std::ostream& f);

namespace auxi {
//! write a dictionnary of things
//! \param dictThings : map of things
//! \param f : output flux
template<MeshMethod meshMethod, class DICT_THING>
void writeDict(string nameObject, const DICT_THING& dictThings, std::ostream& f);
//! write a vector of things
//! \param vecThings : vector of things
//! \param f : output flux
template<MeshMethod meshMethod, class VEC_THING>
void writeVect(string nameObject, VEC_THING vecThings, std::ostream& f);

//! write a point
//! \param object : 2 or 3-D Point
//! \param f : output flux
template<unsigned short DIM, MeshMethod meshMethod>
void write_point(const geoObjects::GeoPoint<DIM>& object, std::ostream& f);
//! write a periodic surface
//! \param object : periodic surface
//! \param f : output flux
template<unsigned short DIM, MeshMethod meshMethod>
void write_perSurface(const geoObjects::PerSurface<DIM>& object, std::ostream& f);
//! write a simple entity (uniquely determined by its name and its leaves)
//! \param nameObject : gmsh name of the object
//! \param object : geometrical entity
//! \param f : output flux
template<MeshMethod meshMethod, class OBJ>
void write_simple(const OBJ& object, std::string nameObject, std::ostream& f);
//! write the sphere in a gmsh format
//! \param f : output flux
template<unsigned short DIM, MeshMethod meshMethod>
void write_sphere(Sphere<DIM> sphere, std::ostream& f, vector<Identifier> removedVolumes = {});
//! write the GeoPerStrucutre<DIM> in a gmsh format
//! \param f : output flux
template<unsigned short DIM, MeshMethod meshMethod>
void write_geoPerStructure(const meshStructure::VoroMesh_Periodic<DIM>& geoPerStructure, std::ostream& f);

//! write the preamble of a mesh description a gmsh format (.geo)
//! \param f : output flux
template<MeshMethod meshMethod>
void writePreamble(std::ostream& f);
//! write the end of a mesh description in a gmsh format (.geo)
//! \param f : output flux
//! \param meshOrder : order of the mesh
//! \param binaryOutput : is the output in binary format?
void writeEnd(std::ostream& f, vector<string> name_outputs = {}, double meshSize = MESH_SIZE, size_t meshOrder = 2, bool binaryOutput = false);
}  // namespace auxi

}  // namespace gmsh_writer
}  // namespace mesh
}  // namespace merope


#include "../Mesh/GmshWriter.ixx"


