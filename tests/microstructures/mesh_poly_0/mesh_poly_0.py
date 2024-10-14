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


L = [1, 1, 1]
nbSpheres = 16 
distMin = 0.05
randomSeed = 0

theSpheres = sac_de_billes.fillMaxRSA_3D(sac_de_billes.Tore, L, nbSpheres, randomSeed, distMin)

for sphere in theSpheres:
    sphere.phase = 2

sphInc = merope.LaguerreTess_3D(L,theSpheres)

mi = merope.MultiInclusions_3D()
mi.setInclusions(sphInc)
mi.addLayer(mi.getAllIdentifiers(), 1, 0.05)


meshGenerator = merope.mesh.MeshGenerator()
meshGenerator.setMeshOrder(2)
meshGenerator.setMeshSize(0.025)
meshGenerator.setMultiInclusions(mi)
meshGenerator.set_nameOutput(["poly.vtk"])
meshGenerator.write("mesh_poly.geo")
