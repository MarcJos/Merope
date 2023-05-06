#!/usr/bin/python
# -*- coding: utf-8 -*- 
# 
# benchmark.py
# Léo Moutin
# 27/10/2021
#
# Compare performances of Neper and mérope
#
# Copyright : see License.txt


import os
import merope
import sac_de_billes
from math import *

import time as time


class resultTest:
    def __init__(self):
        timeList = []
        nbSpheres = 0

def testVoxel(L0, nVox0, nSpheres):
    # Medium dimensions
    L = [L0, L0, L0]
    nVox = [nVox0, nVox0, nVox0]
    N = nSpheres # number of cells
    File_Sphere = "Sphere.txt"

    timeList = []

    ##
    # Generation de la tesselation
    ################################


   
    ##
    # Microstructure Mérope 3D classique
    ################################
    tic0 = time.time()
    r_cel = L0/(4*pi*N/(3*0.38118))**(1./3.)
    spheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, 0, [[r_cel, 1.0]], [1], 0.01) # (sac_de_billes.TypeAlgo, Shape : Tore = periodic cube, seed, RPhi, TabPhase, minDist)
    
    ### building the polyCrystal
    typeCrystal = merope.TypeCrystal.Voronoi       # type of MultiInclusions
    paramMInclusions = merope.SimpleMultiInclusions_3D()
    paramMInclusions.setLength(L)                   # dimensions of the box
    paramMInclusions.setSpheres(spheres)
    paramMInclusions.setTypeCrystal(typeCrystal)
    multiInclusionsCell = paramMInclusions.build()   # build the multiInclusions_3D
    multiInclusionsCell.changePhase(multiInclusionsCell.getAllIdentifiers(),[0 for i in multiInclusionsCell.getAllIdentifiers()])
    multiInclusionsCell.addLayer(multiInclusionsCell.getAllIdentifiers(), 1, 1.0)

    ### Crystal voxellation
    grid = merope.Voxellation_3D(multiInclusionsCell)
    grid.setPureCoeffs([i for i in multiInclusionsCell.getAllIdentifiers()])
    grid.proceed(nVox)
    tic = time.time() - tic0
    timeList.append(tic)
    print("t0 = ", tic)


##
    # Microstructure Mérope 3D simi 2D
    ################################
    tic0 = time.time()
    r_cel = L0/(4*pi*N/(3*0.38118))**(1./3.)
    spheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, 0, [[r_cel, 1.0]], [1], 0.01) # (sac_de_billes.TypeAlgo, Shape : Tore = periodic cube, seed, RPhi, TabPhase, minDist)
    
    ### building the polyCrystal
    typeCrystal = merope.TypeCrystal.Voronoi       # type of MultiInclusions
    paramMInclusions = merope.SimpleMultiInclusions_3D()
    paramMInclusions.setLength(L)                   # dimensions of the box
    paramMInclusions.setSpheres(spheres)
    paramMInclusions.setTypeCrystal(typeCrystal)
    multiInclusionsCell = paramMInclusions.build()   # build the multiInclusions_3D
    multiInclusionsCell.changePhase(multiInclusionsCell.getAllIdentifiers(),[0 for i in multiInclusionsCell.getAllIdentifiers()])
    multiInclusionsCell.addLayer(multiInclusionsCell.getAllIdentifiers(), 1, 1.0)

    ### Crystal voxellation
    grid = merope.Voxellation_3D(multiInclusionsCell)
    grid.setPureCoeffs([i for i in multiInclusionsCell.getAllIdentifiers()])
    nVoxMin = [0, 0, 0]
    nVoxMax = [nVox0, nVox0, 1] # Z slice
    #nVoxMax = [1, nVox0, nVox0] # X slice
    grid.proceed(nVox, nVoxMin, nVoxMax)
    tic = time.time() - tic0
    print("t1 = ", tic)
    timeList.append(tic)


    result = resultTest()
    result.timeList = timeList
    result.nbSpheres = len(spheres)

    return result

