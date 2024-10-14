# -*- coding:utf8 -*-
#
# Tests Multi_Inclusions and Voxellation
# Author: M. Josien
# Date: 10/06/2021
#
# Copyright : see License.txt
#
# Define a recursive structure

import sac_de_billes
import merope

import time

L = [10,10,2]
nWidth = 128
nHeight = 32 


# define structure 1
#--------------------------------------------------------------
## Generate the center of hexagones

i_max = 10
j_max = 10
k_max = 2

spheres = []
for i in range(0,i_max):
    for j in range(0,j_max):
        for k in range(0,k_max):
            x = i + 0.0001
            y = j + 0.5 * (i%2) + 0.0001
            z = k + 0.0001
            spheres.append(sac_de_billes.Sphere_3D([x, y, z], 0, 0))

## building the polyCrystal

polyCrystal = merope.LaguerreTess_3D(L,spheres) 
mIncl = merope.MultiInclusions_3D()
mIncl.setInclusions(polyCrystal)
mIncl.addLayer(mIncl.getAllIdentifiers(), 1000, 0.05)

## define the structure
structure1 = merope.Structure_3D(mIncl)

### define structure 2
#--------------------------------------------------------------
# Define the spheres for the structure #
seedRandom = 0                  # seed for the random generator
desiredRPhi = [[0.5, 1]]  # list of desired[radius,volumeFraction)] for the spheres
tabPhases = [1]           # phases relative to the radii
minDist = 0.                  # minimal distances between the spheres
spheres2 = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seedRandom, desiredRPhi, tabPhases, minDist)
# Define the phase of each inclusions (here, all different)
for i, sphere in enumerate(spheres2):
    sphere.phase = i

polyCrystal2 = merope.LaguerreTess_3D(L, spheres2)
mIncl2 = merope.MultiInclusions_3D()
mIncl2.setInclusions(polyCrystal2)

## define the structure
structure2 = merope.Structure_3D(mIncl2)

### define Mask
seedRandom = 0                  # seed for the random generator
desiredRPhi = [[1, 0.2]]  # list of desired[radius,volumeFraction)] for the spheres
tabPhases = [1]           # phases relative to the radii
minDist = 0.01                  # minimal distances between the spheres
spheresMask = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seedRandom, desiredRPhi, tabPhases, minDist)

sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength(L)
sphIncl.setSpheres(spheresMask)
mInclMask = merope.MultiInclusions_3D()
mInclMask.setInclusions(sphIncl)

## define the mask
structureMask = merope.Structure_3D(mInclMask)

# Combine structures
structureFinale = merope.Structure_3D(structure1, structure2, structureMask)

# Print

def printStructure(structure, name_file_vtk):
    gridParameters = merope.vox.create_grid_parameters_N_L_3D([nWidth,nWidth,nHeight], L)
    grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)
    
    my_printer = merope.vox.vtk_printer_3D()
    my_printer.printVTK_removeUnusedPhase(grid, name_file_vtk, "Coeffs.txt")

printStructure(structure1, "Zone_struct_1.vtk")

printStructure(structure2, "Zone_struct_2.vtk")

printStructure(structureMask, "Zone_struct_mask.vtk")

printStructure(structureFinale, "Zone_struct_finale.vtk")

