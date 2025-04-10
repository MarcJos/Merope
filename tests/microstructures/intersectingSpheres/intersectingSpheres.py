# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/07/2024
#
# Copyright : see License.txt
#
# Create a configuration of penny-shape inclusions
# All are inside non-intersecting spheres
# All are oriented in the same direction

import sac_de_billes
import merope
import sac_de_billes
import time
import numpy as np

sv = 0  # Random number generations seed for separating spheres positioning
N = 150  # Nb spheres
l3D = 1  # RVE dimensions
DIM = 3
L = [l3D, l3D, l3D]
n3D = 64
nbNodes = [n3D, n3D, n3D]

theSpheres = []
theSpheres.append(sac_de_billes.Sphere_3D([0., 0.25, 0.5], 0.25, 1))
theSpheres.append(sac_de_billes.Sphere_3D([0., 0.5, 0.5], 0.25, 2))
theSpheres.append(sac_de_billes.Sphere_3D([0., 0.4, 0.65], 0.25, 3))

e = 0.125  # width of the layer
si = merope.SphereInclusions_3D()
si.setLength(L)
si.setSpheres(theSpheres)

mi = merope.MultiInclusions_3D()
mi.setInclusions(si)
mi.setMatrixPhase(4)
structure = merope.Structure_3D(mi)

gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbNodes, L)
map_phases = {}
for i in range(6):
    for j in range(6):
        map_phases[(i, j)] = 5
gridRepr = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.PolyGeom, map_phases)
gridRepr.convert_to_Iso_format()
gridRepr.apply_coefficients([0., 1., 2., 3., 4., 5.])

print("###########")
print(gridRepr)
print("###########")
gridRepr.apply_homogRule(merope.HomogenizationRule.Largest)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(gridRepr, "myFile.vtk")