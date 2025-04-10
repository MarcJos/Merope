//! Copyright : see license.txt
//!
//! \brief
//
#pragma once


#include "../../../GenericMerope/StdHeaders.hxx"

#include "../Mesh/GeoObjects.hxx"
#include "../Mesh/GmshWriter.hxx"
#include "../Mesh/MeshStructure.hxx"
#include "../Voronoi/VoroToMeshGraph.hxx"


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
//! meshes a microstructure made of Laguerre cells (buggy at the beginning)
//! test for the special functions    setAdimMergeDistance0 // setAdimMergeDistance1(1.e-2);
void test5();
//! meshes a microstructure made of spherical inclusions with holes inside
void test6(int nb_sph = 16);

}  // namespace testMesh
}  // namespace merope




