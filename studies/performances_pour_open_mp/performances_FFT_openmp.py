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
from math import *
import time

output_fileName = "Perf_FFT.txt"

def execution(nb_proc):
    print("###############")
    print(nb_proc)
    print("###############")

    merope.setNbOfThreads(nb_proc)
    tic0 = time.time()
    
    n3D = 64
    L = [10, 10, 10]
    
    
    covariance = lambda x : exp(-x[0]**2 - x[1]**2 - x[2]**2)
    nonLin = lambda x : exp(-x**2)
    gaussianne = merope.gaussianField.GaussianField_3D(covariance, nonLin)
    gaussianne.seed = 1
    
    structure = merope.Structure_3D(gaussianne, L)
    vox = merope.Voxellation_3D(structure)
    vox.proceed([n3D, n3D, n3D])
    
    
        
    final_time = time.time() - tic0
   
    print("##############")
    to_be_printed = str(nb_proc) + " " + str(final_time) + "\n"
    print(to_be_printed)
    print("##############")
    with open(output_fileName, "a") as fic:
        fic.write(to_be_printed)

execution(1)
execution(2)
execution(4)
execution(8)
execution(16)
execution(32)

