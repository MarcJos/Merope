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

n3D = 64
L = [10, 10, 10]


covariance = lambda x : exp(-x[0]**2 - x[1]**2 - x[2]**2)
nonLin = lambda x : exp(-x**2)
gaussianne = merope.gaussianField.GaussianField_3D(covariance, nonLin)
gaussianne.seed = 1

cartesianField = merope.CartesianField_3D(gaussianne, L)
structure = merope.FieldStructure_3D(cartesianField)

vox = merope.Voxellation_3D(structure)
vox.proceed([n3D, n3D, n3D])

vox.printFile("GaussianField.vtk", "Coeffs.txt")

### get the Gaussian Field as an numpy.array
gaussianField = vox.getField()
