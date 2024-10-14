# -*- coding:utf8 -*-
#
# 2D inclusions with a layer
# Author: M. Josien
# Date: 10/06/2021

import sac_de_billes
import merope

# Define the dimensions

L       = [20, 20]
nVox    = [256, 256]

# Throw the seeds
sphIncl = merope.SphereInclusions_2D()
sphIncl.setLength(L)
sphIncl.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[1,1]], [1])
sphIncl.printVER("Seeds.txt")

# Define the polyCrystal
polyCrystal = merope.LaguerreTess_2D(L, sphIncl.getSpheres())

# Define the MultiInclusions

mIncl = merope.MultiInclusions_2D()
mIncl.setInclusions(polyCrystal)

# Define the phases

## Some inclusions have 3 layers, the others have 4 layers
allIdentifiers = mIncl.getAllIdentifiers()
NbInc = len(allIdentifiers)
ThreeLayersIncId    = []
FourLayersIncId     = []
for i in range(0, NbInc):
    if i % 2 == 0:
        ThreeLayersIncId.append(i)
    else:
        FourLayersIncId.append(i)

newPhaseNumber = NbInc
for i in range(0,3):
    mIncl.addLayer(ThreeLayersIncId, newPhaseNumber, 0.15)
    newPhaseNumber += 1

for i in range(0,4):
    mIncl.addLayer(FourLayersIncId, newPhaseNumber, 0.1)
    newPhaseNumber += 1


# Build the grid
structure = merope.Structure_2D(mIncl)

gridParameters = merope.vox.create_grid_parameters_N_L_2D(nVox, L)
gridRepresentation = merope.vox.GridRepresentation_2D(structure, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_2D()
my_printer.printVTK_removeUnusedPhase(gridRepresentation, "Zone.vtk","Coeffs.txt")

