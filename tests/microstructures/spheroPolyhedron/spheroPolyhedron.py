# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 09/09/2022
#
# Copyright : see License.txt
#
# SpheroPolyhedron

import sac_de_billes
import merope
from math import *


L = [4., 4., 4.]
nbVox = [64, 64, 64]

#spheroPolyhedron = merope.microInclusion.SpheroPolyhedronFactory_3D.fromVertices(int(1), [[0., 0., 0.], [1., 0, 0.], [0, 1., 0], [0., 0., 1.]], [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]], float(0.5))
A1 = [1.5, 1.5, 1.5]
A2 = [2.5, 1.5, 1.5]
A3 = [1.5, 2.5, 1.5]
A4 = [1.5, 1.5, 2.5]
faceVert =[[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]
minkowskiRadius = 0.5
spheroPolyhedron = merope.microInclusion.SpheroPolyhedronFactory_3D().fromVertices(int(1), [A1, A2, A3, A4], faceVert, minkowskiRadius)

multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions([spheroPolyhedron], L)
structure = merope.Structure_3D(multiInclusions)

gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbVox, L)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)
grid.apply_coefficients([1, 2])

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Zone.vtk","Coeffs.txt")
