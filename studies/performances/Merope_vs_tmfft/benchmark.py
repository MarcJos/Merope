#!/usr/bin/python
# -*- coding: utf-8 -*- 
# 
# benchmark.py
# Marc Josien
# 01/10/2021
#
# Compare performances of tmfft and mérope
#
# Copyright : see License.txt


import tmfft as fft
import merope
import sac_de_billes

import time as time

class resultTest:
    def __init__(self):
        timeList = []
        nbSpheres = 0

def testVoxel(nbThreads, L0, nVox0, use_tmfft = True):
    fft.setNbOfThreads(nbThreads)
    # Medium dimensions
    L = [L0, L0, L0]
    nVox = [nVox0, nVox0, nVox0]

    File_Sphere = "Sphere.txt"

    timeList = []

    ##
    # Generation sac_de_billes.algo
    ################################
    seed = 0
    desiredRPhi = [[10, 0.2], [5, 0.1]]
    tabPhases = [1,2]
    mindist = 0
    theSpheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seed, desiredRPhi, tabPhases, mindist)
    
    sphereManip = sac_de_billes.SphereManipulator_3D(theSpheres, L)
    sphereManip.printVER(File_Sphere)
    
    ##
    # Microstructure Mérope
    ################################
    tic0 = time.time()
    mIncl =  merope.SimpleMultiInclusions_3D()
    mIncl.setTypeCrystal(merope.TypeCrystal.Voronoi)
    mIncl.setSpheres(theSpheres)
    mIncl.setLength(L)
    
    voxMerope = merope.Voxellation_3D(mIncl.build())
    voxMerope.proceed(nVox)
    
    tic = time.time() -tic0
    timeList.append(tic)
    
    ##
    # Microstructure tmfft
    ################################
    if use_tmfft: 
        tic0 = time.time()
        voro = fft.Voronoi(File_Sphere)
        grid = fft.Grid(voro, nVox[0], nVox[1], nVox[2])
        tic = time.time() - tic0
        timeList.append(tic)
    else:
        timeList.append(0)
    print(timeList)

    result = resultTest()
    result.timeList = timeList
    result.nbSpheres = len(theSpheres)

    return result
