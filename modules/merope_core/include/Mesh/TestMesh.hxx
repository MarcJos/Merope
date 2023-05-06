//! Copyright : see license.txt
//!
//! \brief 
//
#ifndef TESTMESH_HXX_
#define TESTMESH_HXX_


#include "../../../AlgoPacking/src/StdHeaders.hxx"

#include "../Mesh/GeoObjects.hxx"
#include "../Mesh/GmshWriter.hxx"
#include "../Mesh/MeshStructure.hxx"
#include "../Voronoi/VoroToMeshGraph.hxx"

#include "../MeropeNamespace.hxx"


namespace merope {
namespace testMesh {

//!
void test1();
//! tests the interface with voro++ and gmsh
void test2();
//! meshes a microstructure made of spherical inclusions
void test3();
//! meshes a microstructure made of Laguerre cells
void test4();

} // namespace testMesh
} // namespace merope



#endif /* TESTMESH_HXX_ */
