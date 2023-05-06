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

L = [10, 10, 10]
nbVox = [64, 64, 64]

fonction = lambda x : cos(x[0]) * cos(x[1]) * cos(x[2])

scalarField = merope.ScalarField_3D(fonction)
cartesianField = merope.CartesianField_3D(scalarField, L)
structure = merope.FieldStructure_3D(cartesianField)
voxellation = merope.Voxellation_3D(structure)
voxellation.proceed([100, 100, 100])
voxellation.printFile("Field.vtk", "Coeffs.txt")
