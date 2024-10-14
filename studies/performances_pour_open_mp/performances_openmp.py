# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/09/2021
#
# Copyright : see License.txt
#
# Performance tests

import os
import time

import sac_de_billes
import merope

import time


def execution(nb_proc):
    print("###############")
    print(nb_proc)
    print("###############")

    merope.setNbOfThreads(nb_proc)
    tic0 = time.time()
    
    L = [2, 2, 2]
    nVox = [256, 256, 256]
    
    mIncl = merope.SphereInclusions_3D()
    mIncl.setLength(L)
    mIncl.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[0.25,1]], [1])
    mIncl.printVER("Seeds.txt")
    
    polyCrystal = merope.LaguerreTess_3D(L, mIncl.getSpheres())
    mIncl = merope.MultiInclusions_3D()
    mIncl.setInclusions(polyCrystal)
    mIncl.changePhase(mIncl.getAllIdentifiers(),[0 for i in mIncl.getAllIdentifiers()])
    
    mIncl.addLayer(mIncl.getAllIdentifiers(), 1, 0.05)
    mIncl.addLayer(range(0,20), 2 , 0.1)
    mIncl.addLayer(range(10,30), 3, 0.15)
    
    
    print(time.time() - tic0)
    tic0 = time.time()
    
    structure = merope.Structure_3D(mIncl)
    gridParameters = merope.vox.create_grid_parameters_N_L_3D(nVox, L)
    grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Average)
    coefficients = [1. * i for i in mIncl.getAllPhases()]
    grid.apply_homogRule(merope.HomogenizationRule.Voigt, coefficients)
        
    final_time = time.time() - tic0
   
    print("##############")
    to_be_printed = str(nb_proc) + " " + str(final_time) + "\n"
    print(to_be_printed)
    print("##############")
    with open("Perf.txt", "a") as fic:
        fic.write(to_be_printed)

execution(1)
execution(2)
execution(4)
execution(8)
execution(16)
execution(32)

