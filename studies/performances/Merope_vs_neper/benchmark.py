#!/usr/bin/python
# -*- coding: utf-8 -*- 
# 
# benchmark.py
# Marc Josien, adapted by Leo Moutin
# 01/10/2021
#
# Compare performances of Neper and mérope
#
# Copyright : see License.txt


import os
import tmfft as fft
import merope
import sac_de_billes

from reader_tess import *

import time as time

class resultTest:
    def __init__(self):
        timeList = []
        nbSpheres = 0

def testVoxel(L0, nVox0, nSpheres):
    # Medium dimensions
    L = [L0, L0, L0]
    nVox = [nVox0, nVox0, nVox0]
    n = nSpheres # number of cells
    File_Sphere = "Sphere.txt"

    timeList = []

    ##
    # Generation de la tesselation
    ################################

    nameFileNeper = "Neper" + str(L0) + "_" + str(nVox0) + "_"  + str(nSpheres) + ".vtk"
    nameFileMerope = "Merope" + str(L0) + "_" + str(nVox0) + "_"  + str(nSpheres) + ".vtk"


    os.system("neper -T -n " + str(n) +" -domain 'cube("+ str(L0) +","+ str(L0) +","+ str(L0) +")' -morpho 'voronoi' -periodicity all -id 1")
    os.system("mv n"+ str(n) +"-id1.tess Tess_benchmark.tess")
    ##
    # Microstructure Neper
    ################################
    
    tic0 = time.time()
    command = "neper -T -loadtess Tess_benchmark.tess -tesrsize '"+ str(nVox0) +"' -format 'vtk'" 
    print(command)
    os.system(command)

    os.system("mv Tess_benchmark.vtk " + nameFileNeper)
    tic = time.time() - tic0
    timeList.append(tic)
    ##
    # Microstructure Mérope
    ################################
    reader_tess('Tess_benchmark')

    sphereInc = merope.SphereInclusions_3D()
    sphereInc.fromFile("FichierSeeds.txt")
    theSpheres = sphereInc.getSpheres()

	
    ### building the polyCrystal
    tic0 = time.time()
    typeCrystal = merope.TypeCrystal.Voronoi       # type of MultiInclusions
    paramMInclusions = merope.SimpleMultiInclusions_3D()
    paramMInclusions.setLength(L)  
    paramMInclusions.setSpheres(theSpheres)
    paramMInclusions.setTypeCrystal(typeCrystal)
    multiInclusionsCell = paramMInclusions.build()   # build the multiInclusions_3D

    allIdentifiers = multiInclusionsCell.getAllIdentifiers()
    multiInclusionsCell.changePhase(allIdentifiers,[i for i in allIdentifiers])

    ### Crystal voxellation
    grid = merope.Voxellation_3D(multiInclusionsCell)
    grid.setPureCoeffs([i for i in allIdentifiers])
    grid.setVoxelRule(merope.VoxelRule.Center)
    grid.proceed(nVox)
    grid.printFile(nameFileMerope, "Coeffs.txt")
    tic = time.time() - tic0
    timeList.append(tic)
    
    result = resultTest()
    result.timeList = timeList
    result.nbSpheres = len(theSpheres)

    return result
