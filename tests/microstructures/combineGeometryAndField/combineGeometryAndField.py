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


### Define of the micro-structure
def timer_loc(message, t0):
    print(message + " : " + str(time.time() - t0))
    return time.time()

t0 = time.time()

L = [10, 10, 10]
nbVox = [32, 32, 32]

sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength(L)
sphIncl.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[4,0.2],[2,0.2]], [1,2])

multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(sphIncl)

structure = merope.Structure_3D(multiInclusions)

gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbVox, L)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Average)
grid.apply_homogRule(merope.HomogenizationRule.Voigt, [-1,2,3])

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Inclusions.vtk","Coeffs.txt")

t0 = timer_loc("Etape 1", t0)


### Define 3 fields
##
geometryField = grid.get_PureRealField()
cartesianField0 = merope.vox.CartesianField_3D(geometryField, L)
fieldStructure0 = merope.FieldStructure_3D(cartesianField0)
t0 = timer_loc("Etape 2", t0)
#
gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbVox, L)
grid = merope.vox.GridRepresentation_3D(fieldStructure0, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid,"Field0.vtk", "Coeffs0.txt")

t0 = timer_loc("Etape 3", t0)
##
field_function = lambda x : exp(-0.1 * (x[0]**2 + x[1]**2 + x[2]**2))
#scalarField  = merope.ScalarField_3D(field_function) ## deprecated
field_function_pointer = function_py2c(field_function, dim_from=3)
scalarField = merope.ScalarField_3D(field_function_pointer.get_funcPointer())
cartesianField1 = merope.vox.CartesianField_3D(scalarField, L)
fieldStructure1 = merope.FieldStructure_3D(cartesianField1)
t0 = timer_loc("Etape 4", t0)
#
grid = merope.vox.GridRepresentation_3D(fieldStructure1, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Field1.vtk", "Coeffs1.txt")
t0 = timer_loc("Etape 5", t0)
##
covariance = lambda x : exp(-x[0]**2 - x[1]**2 - x[2]**2)
nonLin = lambda x : exp(-x**2)
covariance_pointer = function_py2c(covariance, dim_from=3)
nonLin_pointer = function_py2c(nonLin, dim_from=1)
gaussianne = merope.gaussianField.GaussianField_3D(covariance_pointer.get_funcPointer(), 
        nonLin_pointer.get_funcPointer())
gaussianne.seed = 1
cartesianField2 = merope.vox.CartesianField_3D(gaussianne, L)
fieldStructure2 = merope.FieldStructure_3D(cartesianField2)
t0 = timer_loc("Etape 6", t0)
#
grid = merope.vox.GridRepresentation_3D(fieldStructure2, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Field2.vtk", "Coeffs2.txt")
t0 = timer_loc("Etape 7", t0)
##


### Define fourth field
combine_fieldStructure = merope.FieldStructure_3D(fieldStructure1, fieldStructure2, fieldStructure0)
grid = merope.vox.GridRepresentation_3D(combine_fieldStructure, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "CombineField.vtk", "Coeff3.txt")
t0 = timer_loc("Etape 8", t0)
