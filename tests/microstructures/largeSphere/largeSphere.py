# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/07/2024
#
# Copyright : see License.txt
#
# The sphere is larger than the grid

import sac_de_billes
import merope
import sac_de_billes
import time
import numpy as np


n3D = 32
nbNodes = [n3D, n3D, n3D]

l3D = 1  # RVE dimensions
L = [l3D, l3D, l3D]

theSpheres = []
theSpheres.append(sac_de_billes.Sphere_3D([0., 0., 0.], 1, 1))
theSpheres.append(sac_de_billes.Sphere_3D([0., 0., 0.], 0.5, 2))
#theSpheres.append(sac_de_billes.Sphere_3D([0., 0.4, 0.65], 0.25, 3))

n_spheres = len(theSpheres)

si = merope.SphereInclusions_3D()
si.setLength(L)
si.setSpheres(theSpheres)

mi = merope.MultiInclusions_3D()
mi.setInclusions(si)
mi.setMatrixPhase(0)
structure = merope.Structure_3D(mi)

gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbNodes, L)

gridRepr = merope.vox.GridRepresentation_3D(structure, 
                                            gridParameters, 
                                            merope.vox.VoxelRule.Average)

print("Built the grid")
lambda_sol = 5 
lambda_con = 10
gridRepr.convert_to_Iso_format()
gridRepr.apply_coefficients([0., lambda_sol, lambda_sol, lambda_con])

print("###########")
print(gridRepr)
print("###########")

gridRepr.apply_homogRule(merope.HomogenizationRule.Voigt)
my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(gridRepr, "myFile.vtk")

