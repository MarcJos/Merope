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
grid = merope.Voxellation_3D(multiInclusions)
grid.setPureCoeffs([-1,2,3])
grid.setHomogRule(merope.HomogenizationRule.Voigt)
grid.setVoxelRule(merope.VoxelRule.Average)
grid.proceed(nbVox)

grid.printFile("Inclusions.vtk","Coeffs.txt")
t0 = timer_loc("Etape 1", t0)


### Define 3 fields
##
geometryField = grid.getField()
cartesianField0 = merope.CartesianField_3D(geometryField, L)
fieldStructure0 = merope.FieldStructure_3D(cartesianField0)
t0 = timer_loc("Etape 2", t0)
#
voxellation0 = merope.Voxellation_3D(fieldStructure0)
voxellation0.proceed(nbVox)
voxellation0.printFile("Field0.vtk", "Coeffs0.txt")
t0 = timer_loc("Etape 3", t0)
##
scalarField  = merope.ScalarField_3D(lambda x : exp(-0.1 * (x[0]**2 + x[1]**2 + x[2]**2)))
cartesianField1 = merope.CartesianField_3D(scalarField, L)
fieldStructure1 = merope.FieldStructure_3D(cartesianField1)
t0 = timer_loc("Etape 4", t0)
#
voxellation1 = merope.Voxellation_3D(fieldStructure1)
voxellation1.proceed(nbVox)
voxellation1.printFile("Field1.vtk", "Coeffs1.txt")
t0 = timer_loc("Etape 5", t0)
##
covariance = lambda x : exp(-x[0]**2 - x[1]**2 - x[2]**2)
nonLin = lambda x : exp(-x**2)
gaussianne = merope.gaussianField.GaussianField_3D(covariance, nonLin)
gaussianne.seed = 1
cartesianField2 = merope.CartesianField_3D(gaussianne, L)
fieldStructure2 = merope.FieldStructure_3D(cartesianField2)
t0 = timer_loc("Etape 6", t0)
#
voxellation2 = merope.Voxellation_3D(fieldStructure2)
voxellation2.proceed(nbVox)
voxellation2.printFile("Field2.vtk", "Coeffs2.txt")
t0 = timer_loc("Etape 7", t0)
##


### Define fourth field
combine_fieldStructure = merope.FieldStructure_3D(fieldStructure1, fieldStructure2, fieldStructure0)
voxellation3 = merope.Voxellation_3D(combine_fieldStructure)
voxellation3.proceed(nbVox)
voxellation3.printFile("CombineField.vtk", "Coeff3.txt")
t0 = timer_loc("Etape 8", t0)
