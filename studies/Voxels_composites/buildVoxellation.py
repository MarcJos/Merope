# -*- coding:utf8 -*-
#
# Example for the use of composite voxels.
# We want to compute the conductivity of lead spheres coated with gold inside water.
#
# Author: M. Josien
# Date:  2021-12-16
#
#
# Copyright : see License.txt

import os

import sac_de_billes
import merope

import archi_merope as arch
import amitex_fftp.Lancement_AMITEX as amitex
import amitex_fftp.Post_traitement as amitex_out

from parametrization import *

GRID_FILE = "Grid.vtk"

def geometry():
        #--------------------------------------------------------------
    # Define the spheres for the structure #
    typeAlgo = sac_de_billes.TypeAlgo.WP     # method for throwing spheres
    nameShape = sac_de_billes.NameShape.Tore  # periodic cube
    seedRandom = 0                  # seed for the random generator
    desiredRPhi = [[2, 0.65]]
    tabPhases = [1]           # phases relative to the radii
    minDist = 0                  # minimal distances between the spheres
    sphereList = sac_de_billes.throwSpheres_3D(typeAlgo, nameShape, L, seedRandom, desiredRPhi, tabPhases, minDist)
    #----------------------------------------------------------------
    # Define the structure
    sphereInclusions = merope.SphereInclusions_3D()
    sphereInclusions.setLength(L)
    sphereInclusions.setSpheres(sphereList)
    
    multiInclusions = merope.MultiInclusions_3D()
    multiInclusions.setInclusions(sphereInclusions)
    multiInclusions.addLayer(multiInclusions.getAllIdentifiers(), 2, 0.1)
    return multiInclusions


def createGrid(n, voxelRule, homogRule, multiInclusions):
# Voxellate
    nbVox = [n, n, n]
    grid = merope.Voxellation_3D(multiInclusions)
    grid.setVoxelRule(voxelRule)
    grid.setHomogRule(homogRule)
    grid.setPureCoeffs(all_lambdas)
    grid.proceed(nbVox)
    grid.printFile(GRID_FILE, "Coeffs.txt")
    ###
    os.system("rm -rf "     + arch.resultFolder)
    os.system("mkdir "      + arch.resultFolder)
    os.system("mv *txt "    + arch.resultFolder)
    os.system("mv *vtk "    + arch.resultFolder)

def lanceAmitex():
    ### Launch amitex
    os.chdir(arch.resultFolder)
    number_of_processors = 16        ### for parallel computing
    ### Implicitly uses the file "Coeffs.txt" to compute the thermal matrix
    amitex.computeThermalCoeff(GRID_FILE, number_of_processors)
    ### See the homogenized matrix
    for i in range(0,10**6):
        pass
    homogenized_matrix = amitex_out.printThermalCoeff(".")
    with open("thermalCoeff_amitex.txt", "w") as fic:
        for line in homogenized_matrix:
            for coeff in line:
                fic.write(str(coeff))
                fic.write(" ")
            fic.write("\n")
    os.chdir("../")
    return homogenized_matrix

def scalarCoeff(matrix):
    result = 0
    for i in range(0,3):
        result += matrix[i][i]
    return result /3

def wholeProcedure(n, voxelRule, homogRule, multiInclusions):
    createGrid(n, voxelRule, homogRule, multiInclusions)
    return scalarCoeff(lanceAmitex())
