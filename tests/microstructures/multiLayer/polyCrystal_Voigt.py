# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 10/06/2021
# 
# Copyright : see License.txt
#
# Build a polycrystal with layer

import sac_de_billes
import merope

import time

tic0 = time.time()

L = [2, 2, 2]
nVox = [100, 100, 100]

mIncl = merope.SphereInclusions_3D()
mIncl.setLength(L)
mIncl.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[0.25,1]], [1])
mIncl.printVER("Seeds.txt")

polyCrystal = merope.LaguerreTess_3D(L, mIncl.getSpheres())
mIncl = merope.MultiInclusions_3D()
mIncl.setInclusions(polyCrystal)
mIncl.changePhase(mIncl.getAllIdentifiers(),[0 for i in mIncl.getAllIdentifiers()])

mIncl.addLayer(mIncl.getAllIdentifiers(), 1, 0.05)
mIncl.addLayer(range(0,20), 2 , 0.1)
mIncl.addLayer(range(10,30), 3, 0.15)


print(time.time() - tic0)
tic0 = time.time()

grid = merope.Voxellation_3D(mIncl)
grid.setVoxelRule(merope.VoxelRule.Average)
grid.setHomogRule(merope.HomogenizationRule.Voigt)
grid.proceed(nVox)
grid.printFile("Zone_Voigt.vtk","Coeffs_Voigt.txt")

print(time.time() - tic0)
