//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef MESH_GMSHWRITER_HXX_
#define MESH_GMSHWRITER_HXX_

#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Mesh/GeoObjects.hxx"
#include "../Mesh/MeshStructure.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {

namespace mesh {
namespace gmsh_writer {

inline string mesh_size = "meshSize";
inline constexpr double MESH_SIZE = 0.05;
inline constexpr size_t NUMBER_MESHCOMPONENT_PER_SPHERE = 15;

//! write a geometrical entity in a gmsh format (.geo)
//! \param object : geometrical entity
//! \param f : output flux
template<class OBJ>
void write(const OBJ& object, std::ostream& f);

namespace auxi {
//! write a dictionnary of things
//! \param dictThings : map of things
//! \param f : output flux
template<class DICT_THING>
void writeDict(string nameObject, const DICT_THING& dictThings, std::ostream& f);
//! write a vector of things
//! \param vecThings : vector of things
//! \param f : output flux
template<class VEC_THING>
void writeVect(string nameObject, VEC_THING vecThings, std::ostream& f);

//! write a point
//! \param object : 2 or 3-D Point
//! \param f : output flux
template<unsigned short DIM>
inline void write_point(const geoObjects::GeoPoint<DIM>& object, std::ostream& f);
//! write a periodic surface
//! \param object : periodic surface
//! \param f : output flux
template<unsigned short DIM>
inline void write_perSurface(const geoObjects::PerSurface<DIM>& object, std::ostream& f);
//! write a simple entity (uniquely determined by its name and its leaves)
//! \param nameObject : gmsh name of the object
//! \param object : geometrical entity
//! \param f : output flux
template<class OBJ>
inline void write_simple(const OBJ& object, std::string nameObject, std::ostream& f);
//! write the sphere in a gmsh format
//! \param f : output flux
template<unsigned short DIM>
inline void write_sphere(Sphere<DIM> sphere, std::ostream& f, vector<Identifier> removedVolumes = {});
//! write the GeoPerStrucutre<DIM> in a gmsh format
//! \param f : output flux
template<unsigned short DIM>
inline void write_geoPerStructure(const meshStructure::VoroMesh_Periodic<DIM>& geoPerStructure, std::ostream& f);

//! write the preamble of a mesh description a gmsh format (.geo)
//! \param f : output flux
void writePreamble(std::ostream& f);
//! write the end of a mesh description in a gmsh format (.geo)
//! \param f : output flux
//! \param meshOrder : order of the mesh
//! \param binaryOutput : is the output in binary format?
void writeEnd(std::ostream& f, double meshSize = MESH_SIZE, size_t meshOrder = 2, bool binaryOutput = false);
} // namespace auxi

} // namespace gmsh_writer
} // namespace mesh
} // namespace merope


#include "../Mesh/GmshWriter.ixx"

#endif /* MESH_GMSHWRITER_HXX_ */
