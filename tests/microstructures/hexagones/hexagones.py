# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 10/06/2021
#
# Copyright : see License.txt
#
# Build hexagonal structure


import sac_de_billes
import merope

import time

L = [10,10,2]
n3D =500


### Generate the center of hexagones

i_max = 10
j_max = 10
k_max = 2

spheres = []
for i in range(0,i_max):
    for j in range(0,j_max):
        for k in range(0,k_max):
            x = i + 0.0001
            y = j + 0.5 * (i%2) + 0.0001
            z = k + 0.0001
            spheres.append(sac_de_billes.Sphere_3D([x, y, z], 0, 0))

### building the polyCrystal

polyCrystal = merope.LaguerreTess_3D(L, spheres)
polyCrystal.computeTessels()
mIncl = merope.MultiInclusions_3D()
mIncl.setInclusions(polyCrystal)
mIncl.addLayer(mIncl.getAllIdentifiers(), 1000, 0.05)
structure = merope.Structure_3D(mIncl)

### voxellation
gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D,n3D,10], L)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)
pureCoeffs = [2 for i in range(0,1001)]
pureCoeffs[1000] = 1
grid.apply_coefficients(pureCoeffs)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(grid, "Field.vtk")
my_printer.printVTK_segmented(grid, "Zone.vtk","Coeffs.txt")

### test of structure
structure = merope.Structure_3D(mIncl)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)
pureCoeffs = [2 for i in range(0,1001)]
pureCoeffs[1000] = 1
grid.apply_coefficients(pureCoeffs)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Zone_bis.vtk","Coeffs_bis.txt")

