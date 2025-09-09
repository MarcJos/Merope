# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 04/10/2022
#
# Copyright : see License.txt
#
# combine geometry and fields

import sac_de_billes
import merope
from math import *
import time
from lambda_function_builder import *


tic_0 = time.time()

### Define of the micro-structure
def timer_loc(message, t0):
    print(message + " : " + str(time.time() - t0))
    return time.time()

t0 = time.time()

L = [10, 10, 10]
nbVox = [256, 256, 256]

gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbVox, L)


#### polycrystal
seed = 0
spheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seed, [[2, 1]], [1], 0)

polyCrystal = merope.LaguerreTess_3D(L, spheres)
multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(polyCrystal)
multiInclusions.changePhase(multiInclusions.getAllIdentifiers(),multiInclusions.getAllIdentifiers())

structure = merope.Structure_3D(multiInclusions)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)

"""
my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(grid, "Inclusions.vtk","Coeffs.txt")
"""

t0 = timer_loc("Etape 1", t0)


### Define 3 fields
##
grid.apply_coefficients(multiInclusions.getAllIdentifiers())
geometryField = grid.get_PureRealField()
cartesianField0 = merope.vox.CartesianField_3D(geometryField, L)
fieldStructure0 = merope.FieldStructure_3D(cartesianField0)
t0 = timer_loc("Etape 2", t0)
#
grid = merope.vox.GridRepresentation_3D(fieldStructure0, gridParameters, merope.vox.VoxelRule.Center)

"""
my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid,"Field0.vtk", "Coeffs0.txt")
"""

t0 = timer_loc("Etape 5", t0)
##
covariance = lambda x : exp(- 200 * (x[0]**2 + x[1]**2 + x[2]**2))
nonLin = lambda x : exp(-x**2)
covariance_pointer = function_py2c(covariance, dim_from=3)
nonLin_pointer = function_py2c(nonLin, dim_from=1)
gaussianne = merope.gaussianField.GaussianField_3D(covariance_pointer.get_funcPointer(), 
        nonLin_pointer.get_funcPointer())
gaussianne.seed = 2
cartesianField2 = merope.vox.CartesianField_3D(gaussianne, L)
fieldStructure2 = merope.FieldStructure_3D(cartesianField2)
t0 = timer_loc("Etape 6", t0)
#
grid = merope.vox.GridRepresentation_3D(fieldStructure2, gridParameters, merope.vox.VoxelRule.Center)

"""
my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(grid, "Field2.vtk")
t0 = timer_loc("Etape 7", t0)
"""
##


### Define fourth field

fonction = lambda x : x[0] + x[1]
fonc_pointer = function_py2c(fonction, dim_from=2)
interf_fonc = fonc_pointer.get_funcPointer()

combine_fieldStructure = merope.FieldStructure_3D(fieldStructure0, fieldStructure2, interf_fonc)
grid = merope.vox.GridRepresentation_3D(combine_fieldStructure, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(grid, "CombineField.vtk")
t0 = timer_loc("Etape 8", t0)


print("Total time =", time.time() - tic_0)
