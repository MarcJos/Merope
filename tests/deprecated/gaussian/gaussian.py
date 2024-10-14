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


print("This way of using lambda functions for Gaussian field is deprecated!")

n3D = 64
L = [10, 10, 10]

def printField(field, L, n3D, name_vtk):
  cartesianField = merope.vox.CartesianField_3D(field, L)
  structureG = merope.FieldStructure_3D(cartesianField)

  gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D, n3D, n3D], L)
  grid = merope.vox.GridRepresentation_3D(structureG, gridParameters, merope.vox.VoxelRule.Center)
  
  my_printer = merope.vox.vtk_printer_3D()
  my_printer.printVTK(grid, name_vtk)


def createGaussianField(lambda_multiplication, name, seed):
  covariance = lambda x : exp((-x[0]**2 - x[1]**2 - x[2]**2) * lambda_multiplication)
  nonLin = lambda x : lambda_multiplication *  exp(-x**2)
  gaussianne = merope.gaussianField.GaussianField_3D(covariance, nonLin)
  gaussianne.seed = seed
  
  ### print the desired covariance
  covField = merope.ScalarField_3D(covariance)
  printField(covField, L, n3D, name + "desired_cov.vtk")
  ### print the actual covariance
  seenCovField =  merope.gaussianField.NumericalCovariance_3D(covariance)
  printField(seenCovField, L, n3D, name + "real_cov.vtk")

  cartesianField = merope.vox.CartesianField_3D(gaussianne, L)
  structureG = merope.FieldStructure_3D(cartesianField)
  
  gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D, n3D, n3D], L)
  grid = merope.vox.GridRepresentation_3D(structureG, gridParameters, merope.vox.VoxelRule.Center)
  
  my_printer = merope.vox.vtk_printer_3D()
  my_printer.printVTK_segmented(grid, name + "GaussianField.vtk", name + "Coeffs.txt")
  return structureG
  
structureGaussianField_0 = createGaussianField(1, "", 1)
structureGaussianField_1 = createGaussianField(5, "1_", 0)



### Create a polycrystal
seed = 0
spheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seed, [[2, 1]], [1], 0)
polyCrystal = merope.LaguerreTess_3D(L, spheres)
multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(polyCrystal)
multiInclusions.changePhase(multiInclusions.getAllIdentifiers(), [i%5 for i in multiInclusions.getAllIdentifiers()])
structureMu = merope.Structure_3D(multiInclusions)

gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D, n3D, n3D], L)
grid = merope.vox.GridRepresentation_3D(structureMu, gridParameters, merope.vox.VoxelRule.Average)
grid.apply_homogRule(merope.HomogenizationRule.Voigt)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "poly.vtk", "cpoly.txt")

### Combine structures
#  Get 2nd microstructure as FieldStructure
fieldMu = voxellationMu.get_PureRealField()
cartesianMu = merope.vox.CartesianField_3D(fieldMu, L)
fieldStructureMu = merope.FieldStructure_3D(cartesianMu)

#Combine
fieldStructureTot = merope.FieldStructure_3D(structureGaussianField_0, structureGaussianField_1, fieldStructureMu)

structure = merope.Structure_3D(multiInclusions)

gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D, n3D, n3D], L)
grid = merope.vox.GridRepresentation_3D(fieldStructureTot, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "totalStruct.vtk", "totalStruct.txt")








