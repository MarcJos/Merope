# -*- coding:utf8 -*-
#
# Tests Multi_Inclusions and Voxellation
# Author: M. Josien
# Date: 10/06/2021

import sac_de_billes
import merope

import time

tic0 = time.time()

L = [10, 10, 10]
n3D = 256 

### Add the inclusions

sphIncl2 = merope.SphereInclusions_3D()
sphIncl2.setLength(L)
sphIncl2.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[0.25, 1]], [1])

multiInclusions2 = merope.MultiInclusions_3D()
multiInclusions2.setInclusions(sphIncl2)

### Get the spherical inclusions

grid = merope.Voxellation_3D(multiInclusions2)
grid.setPureCoeffs([0, 1])
grid.setVoxelRule(merope.VoxelRule.Center)
grid.setHomogRule(merope.HomogenizationRule.Largest)
grid.proceed([n3D,n3D,n3D])
grid.printFile("Zone_Inclusions.vtk","Coeffs.txt")

### Laguerre

sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength(L)
sphIncl.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[2.,1.]], [1])
sphIncl.printVER("Seeds.txt")

polyCrystal = merope.LaguerreTess_3D(L, sphIncl.getSpheres())
multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(polyCrystal)
N = len(multiInclusions.getAllIdentifiers())
multiInclusions.addLayer(multiInclusions.getAllIdentifiers(), N, 0.3)
multiInclusions.changePhase(multiInclusions.getAllIdentifiers(), [1 for i in multiInclusions.getAllIdentifiers()])

### Get the polyCrystal

grid = merope.Voxellation_3D(multiInclusions)
grid.setPureCoeffs([i for i in range(0,N+1)])
grid.setVoxelRule(merope.VoxelRule.Center)
grid.setHomogRule(merope.HomogenizationRule.Largest)
grid.proceed([n3D,n3D,n3D])
grid.printFile("Zone_Crystal.vtk","Coeffs.txt")


phase = [N]
dictionnaire = {i : N+1  for i in phase}

structure = merope.Structure_3D(multiInclusions, multiInclusions2, dictionnaire)

###

grid = merope.Voxellation_3D(structure)
grid.setPureCoeffs([i for i in range(0,N+2)])
grid.setVoxelRule(merope.VoxelRule.Center)
grid.setHomogRule(merope.HomogenizationRule.Largest)
grid.proceed([n3D,n3D,n3D])
grid.printFile("Zone.vtk","Coeffs.txt")


