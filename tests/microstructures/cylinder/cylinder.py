# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/07/2024
#
# Copyright : see License.txt
#
# Test cylinder

import sac_de_billes
import merope
import time

L = [10, 10, 10]
axis = merope.geometry.Segment_3D([[5, 9, 5], [3, 12, 5]])
cylinder = merope.geometry.Cylinder(axis, 1)
cylinderInc = merope.microInclusion.CylinderInclusions(cylinder)

cylinderInclusions = merope.CylinderInclusions_3D()
cylinderInclusions.setLength(L)
cylinderInclusions.setInclusions([cylinderInc])

multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(cylinderInclusions)
multiInclusions.setMatrixPhase(1)

### pure case

structure = merope.Structure_3D(multiInclusions)
gridParameters = merope.vox.create_grid_parameters_N_L_3D([100,100,100], L)

grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)
grid.apply_coefficients([0,1,2])

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Zone.vtk","Coeffs.txt")


### composite case

grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Average)
grid.apply_homogRule(merope.HomogenizationRule.Reuss, [1,2])

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Zone_composite.vtk","Coeffs_composite.txt")

