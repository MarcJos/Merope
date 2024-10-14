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


structure = merope.Structure_3D(multiInclusions)

gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D,n3D,n3D], L)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(grid, "Zone.vtk")

print(time.time() - tic0)
tic0 = time.time()

structure = merope.Structure_3D(multiInclusions)

grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Average)
coefficients = [1,2,3, 4, 5, 6, 7, 8, 9]
grid.apply_homogRule(merope.HomogenizationRule.Voigt, coefficients)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Zone_Voigt.vtk","Coeffs_Voigt.txt")

print(time.time() - tic0)
tic0 = time.time()


polyX = merope.SimpleStructure_3D()
polyX.mainInclusions.setTypeCrystal(merope.TypeCrystal.Laguerre)
polyX.mainInclusions.fromFile("Seeds.txt")

# Grid
structure = polyX.build()
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(grid, "Cristal.vtk")

print(time.time() - tic0)
