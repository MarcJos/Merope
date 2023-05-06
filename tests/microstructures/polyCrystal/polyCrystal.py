# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 10/06/2021
#
# Copyright : see License.txt
#
# Build a polycrystal

import sac_de_billes
import merope

import time

tic0 = time.time()

L = [2, 2, 2]
n3D =64


seed = 0
spheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seed, [[0.5, 1]], [1], 0)

### print
mIncl = merope.SphereInclusions_3D()
mIncl.setLength(L)
mIncl.setSpheres(spheres)
mIncl.printVER("Seeds.txt")
###


polyCrystal = merope.LaguerreTess_3D(L, spheres)
multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(polyCrystal)
multiInclusions.changePhase(multiInclusions.getAllIdentifiers(),multiInclusions.getAllIdentifiers())


grid = merope.Voxellation_3D(multiInclusions)
grid.setPureCoeffs([1,2,3, 4, 5, 6, 7, 8, 9])
grid.setVoxelRule(merope.VoxelRule.Center)
grid.setHomogRule(merope.HomogenizationRule.Voigt)
grid.proceed([n3D,n3D,n3D])
grid.printFile("Zone.vtk","Coeffs.txt")

print(time.time() - tic0)
tic0 = time.time()

grid = merope.Voxellation_3D(multiInclusions)
grid.setPureCoeffs([1,2,3, 4, 5, 6, 7, 8, 9])
grid.setVoxelRule(merope.VoxelRule.Average)
grid.setHomogRule(merope.HomogenizationRule.Voigt)
grid.proceed([n3D,n3D,n3D])
grid.printFile("Zone_Voigt.vtk","Coeffs_Voigt.txt")

print(time.time() - tic0)
tic0 = time.time()


polyX = merope.SimpleStructure_3D()
polyX.mainInclusions.setTypeCrystal(merope.TypeCrystal.Laguerre)
polyX.mainInclusions.fromFile("Seeds.txt")

# Grid
g3D = merope.Voxellation_3D(polyX.build())
g3D.proceed([n3D, n3D, n3D])
g3D.printFile("Cristal.vtk", "Coeffs_0.txt")

print(time.time() - tic0)
