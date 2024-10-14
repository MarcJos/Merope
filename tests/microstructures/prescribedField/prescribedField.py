# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 09/09/2022
#
# Copyright : see License.txt
#
# Build Gaussian fields


import sac_de_billes
import merope
from math import *
from lambda_function_builder import *
import time

L = [10, 10, 10]
nbVox = [64, 64, 64]

def pb(nb_proc):
    t = time.time()
    merope.setNbOfThreads(nb_proc) 

    fonction = lambda x : cos(x[0]) * cos(x[1]) * cos(x[2])
    fonc_pointer = function_py2c(fonction, dim_from=3)
    interf_fonc = fonc_pointer.get_funcPointer()
    
    scalarField = merope.ScalarField_3D(interf_fonc)
    cartesianField = merope.vox.CartesianField_3D(scalarField, L)
    structure = merope.FieldStructure_3D(cartesianField)
    
    gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbVox, L)
    grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)

    my_printer = merope.vox.vtk_printer_3D()
    my_printer.printVTK_segmented(grid, "Field.vtk", "Coeffs.txt")

    print("Nb proc : " + str(nb_proc) + " | Time used : " +  str(time.time() - t))



pb(1)
pb(2)
pb(4)
pb(8)
