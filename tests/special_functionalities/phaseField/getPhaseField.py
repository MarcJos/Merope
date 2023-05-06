# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 29/08/2022
#
# Copyright : see License.txt
#
# Get a phase field from a microstructure

# Import Merope & tmfft
import sac_de_billes
import merope
import tmfft as tmfft
import archi_merope as arch

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

nbVox = [100,100,100]

### Play with the phaseField
def printHashtags():
  for i in range(10):
    print("##############")

phaseField = grid.computePhaseGrid(nbVox)

printHashtags()
for i in range(3,8):
  j = 5
  k = 12
  print("The content of the voxel at position in discrete coordinates [i, j ,k] = " + str([i,j,k]) + " is ")
  linear_discrete_coordinate = merope.get_linear_index_3D([i, j, k], nbVox)
  print(phaseField[linear_discrete_coordinate])
  print("which is a list of [phase, volumeFraction]")
printHashtags()

### Print the grid
grid.proceed(nbVox)
grid.printFile("Zone.vtk","Coeffs.txt")

