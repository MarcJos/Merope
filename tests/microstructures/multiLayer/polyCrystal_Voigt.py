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


structure = merope.Structure_3D(mIncl)

gridParameters = merope.vox.create_grid_parameters_N_L_3D(nVox, L)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Average)
coeffs = [i for i in structure.getAllPhases()]
grid.apply_homogRule(merope.HomogenizationRule.Voigt, coeffs)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Zone_Voigt.vtk","Coeffs_Voigt.txt")

print(time.time() - tic0)
