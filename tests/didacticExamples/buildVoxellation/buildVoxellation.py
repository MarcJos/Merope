# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/09/2021
#
# Copyright : see License.txt
#
# Didactic example for :
# - defining the random distribution of sphere for the structure
# - defining a simple structure from a sphere
# - add layers to it
# - voxellize it

# Import Merope & tmfft
import sac_de_billes
import merope
import tmfft as tmfft
import archi_merope as arch

# Dimensions of the box
L = [10., 10., 10.]

#--------------------------------------------------------------
# Define the spheres for the structure #
typeAlgo = sac_de_billes.TypeAlgo.RSA     # method for throwing spheres
nameShape = sac_de_billes.NameShape.Tore  # periodic cube
seedRandom = 0                  # seed for the random generator
desiredRPhi = [[1, 0.1], [0.5, 0.1], [0.25, 0.1]]  # list of desired[radius,volumeFraction)] for the spheres
tabPhases = [0, 1, 2]           # phases relative to the radii
minDist = 0.01                  # minimal distances between the spheres
sphereList = sac_de_billes.throwSpheres_3D(typeAlgo, nameShape, L, seedRandom, desiredRPhi, tabPhases, minDist)
nbInclusions = len(sphereList)
# Define the phase of each inclusions (here, all different)
for i, sphere in enumerate(sphereList):
    sphere.phase = i
#--------------------------------------------------------------

#--------------------------------------------------------------
# Define the main inclusions
typeCrystal = merope.TypeCrystal.Laguerre       # type of MultiInclusions
aspectRatio = [1, 2, 3]                         # aspect ratio of the inclusions
#
paramMInclusions = merope.SimpleMultiInclusions_3D()
paramMInclusions.setLength(L)                   # dimensions of the box
paramMInclusions.setSpheres(sphereList)
paramMInclusions.setTypeCrystal(typeCrystal)
paramMInclusions.setAspRatio(aspectRatio)       # useless for TypeCrystal.Spheres
multiInclusions_3D = paramMInclusions.build()   # build the multiInclusions_3D
#--------------------------------------------------------------

#--------------------------------------------------------------
# Define layers
allIdentifiers = multiInclusions_3D.getAllIdentifiers() # get all the identifiers of the inclusions, in order to modify them
identifiers = allIdentifiers                                            # which idenfied inclusion a layer will be added
newPhase = [nbInclusions + i for i, ident in enumerate(identifiers)]    # phases of the new layers
width = [0.05 for ident in identifiers]                                 # width of the new layers
multiInclusions_3D.addLayer(identifiers, newPhase, width)   # add the new layers
# get Infos:
allPhases = multiInclusions_3D.getAllPhases()               # list all the phases in the structure
#--------------------------------------------------------------

#--------------------------------------------------------------
# Voxellation 1 : print the phases
nbVox = [128, 128, 128]     # nb of voxels in each direction
# Voxellation 1 : print the phases
voxellation_3D = merope.Voxellation_3D(multiInclusions_3D)  # constructor
voxellation_3D.proceed(nbVox)                               # build the grid
voxellation_3D.printFile("Phases.vtk", "Phases.txt")        # print the grid
#--------------------------------------------------------------

#--------------------------------------------------------------
# Voxellation 2 : composite voxels
## Define coefficients
lambdaCoefficients = [1, 2, 3]
pureCoeffs = [lambdaCoefficients[i%3] for i in allPhases]                    # define coefficients (of pure phases)
## Parametrize voxellation
voxelRule = merope.VoxelRule.Average # the volume fraction of each phase is computed in every voxel
homogRule = merope.HomogenizationRule.Voigt # the coefficient are averaged according to that rule
## Build voxellation
voxellation_3D_composite = merope.Voxellation_3D(multiInclusions_3D)   # constructor
voxellation_3D_composite.setVoxelRule(voxelRule)
voxellation_3D_composite.setHomogRule(homogRule)
voxellation_3D_composite.setPureCoeffs(pureCoeffs) # define the coefficients relative to the phases
voxellation_3D_composite.proceed(nbVox)
voxellation_3D_composite.printFile("Composite.vtk", "Coeffs.txt")
#--------------------------------------------------------------

# Move the results to a subdirectory
import os
os.system("mkdir "      + arch.resultFolder)
os.system("mv *txt "    + arch.resultFolder)
os.system("mv *vtk "    + arch.resultFolder)

