# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/09/2021
#
# Copyright : see License.txt
#
# Build a .geo file to be meshed by gmsh
# The geometry is made of polyhedra


import sac_de_billes
import merope


L = [2, 2, 2]
nbSpheres = 16 
distMin = 0.1
randomSeed = 0

theSpheres = sac_de_billes.fillMaxRSA_3D(sac_de_billes.Tore, L, nbSpheres, randomSeed, distMin)

for sphere in theSpheres:
    sphere.phase = 2

sphInc = merope.SphereInclusions_3D()
sphInc.setLength(L)
sphInc.setSpheres(theSpheres)

mi = merope.MultiInclusions_3D()
mi.setInclusions(sphInc)
mi.setMatrixPhase(1)
mi.addLayer(mi.getAllIdentifiers(), 3, 0.1)
mi.addLayer(mi.getAllIdentifiers(), 4, 0.1)


meshGenerator = merope.mesh.MeshGenerator()
meshGenerator.setMeshOrder(2)
meshGenerator.setMeshSize(0.05)
meshGenerator.setMultiInclusions(mi)
meshGenerator.do_not_mesh([2])
meshGenerator.set_nameOutput(["spheres.vtk"])
meshGenerator.write("mesh_spheres.geo")

meshGenerator.write("mesh_spheres_OCC.geo", merope.mesh.MeshMethod.OpenCascade)

