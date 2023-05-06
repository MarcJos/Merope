#!/usr/bin/python
# coding: utf-8
#
# Build overlapping coated inclusions
#

import merope
import sac_de_billes
import time

tic0 = time.time()
n3D = 256       # Number of pixel by dimentions
L = [10,10,10]

R1 = 1  # Radius of the inclusion with coating
R2 = 0.2 # Radius of the central inclusion
layerWidth = R1-R2  # Thickness of the coating

# Definition of the uniformous matrix
sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength(L)
mIncl = merope.MultiInclusions_3D()
mIncl.setInclusions(sphIncl)
mIncl.setMatrixPhase(0)
structure1 = merope.Structure_3D(mIncl)

# Creation of the mask 
seedRandom = 0              # seed for the random generator
desiredRPhi = [[R1, 0.4]]     # list of desired[radius,volumeFraction)] for the spheres
tabPhases = [1]             # phases relative to the radii
minDist = -2*layerWidth     # Minimal distances between the spheres (the coating can percolate but not the inclusion)

spheresMask = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.WP, sac_de_billes.NameShape.Tore, L, seedRandom, desiredRPhi, tabPhases, minDist) 
for sphere in spheresMask:
    sphere.radius = R1
    # For the mask the matrix is phase 0 and the big inclusion (same size as the inclusion with coating) is phase 1 (only phase 1 will be replaced)
    sphere.phase = 1 
nb_sphere = len(spheresMask)
sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength(L)
sphIncl.setSpheres(spheresMask)
mInclMask = merope.MultiInclusions_3D()
mInclMask.setInclusions(sphIncl)
structureMask = merope.Structure_3D(mInclMask)

# Second structure which will be inside the mask
for sphere in spheresMask:
    sphere.radius = R2 # Reduction of the radius (inclusion with no coating)
    sphere.phase = 2 # The central inclusion is phase 2
sphIncl2 = merope.SphereInclusions_3D()
sphIncl2.setLength(L)
sphIncl2.setSpheres(spheresMask)
mIncl2 = merope.MultiInclusions_3D()
mIncl2.setInclusions(sphIncl2)
mIncl2.setMatrixPhase(1) # Necessity to change the phase of the matrix to be different from the matrix of the first structure, if not all the phases will be the same and it results in only 2 phases instead of 3
structure2 = merope.Structure_3D(mIncl2)

# Combination of the structures
structureFinale = merope.Structure_3D(structure1, structure2, structureMask)

# Print
grid = merope.Voxellation_3D(structureFinale)
grid.proceed([n3D,n3D,n3D])
grid.printFile("Zone_struct_finale.vtk","Coeffs_finale.txt")
