# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 10/06/2021
#
# Copyright : see License.txt
#
# Tests Multi_Inclusions and Voxellation

import sac_de_billes
import merope
import time

tic0 = time.time()

### Define of the micro-structure

sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength([10,10,10])
sphIncl.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[0.5,0.2],[0.25,0.2]], [1,2])

### Print outputs

sphIncl.printDump("inclusions.dump")
sphIncl.printPos("inclusions.pos")
sphIncl.printCSV("inclusions.csv")
sphIncl.printFracVol("volume_fractions.txt")

### Build the voxellation

multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(sphIncl)
grid = merope.Voxellation_3D(multiInclusions)
grid.setPureCoeffs([1,2,3])
grid.setHomogRule(merope.HomogenizationRule.Voigt)
grid.setVoxelRule(merope.VoxelRule.Average)
grid.proceed([100,100,100])
grid.printFile("Zone.vtk","Coeffs.txt")

print(time.time() - tic0)
