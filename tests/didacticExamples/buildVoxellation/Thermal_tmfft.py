# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/09/2021
#
# Copyright : see License.txt
#
# Script for launching the fft solver tmfft
#
import os
import archi_merope as arch

import tmfft as fft 
from math import *
import numpy as np


### Launch tmfft
number_of_processors = 16                   ### for parallel computing
voxellation_of_zones = "Composite.vtk"      ### produced by Merope (voxellation)
coeffFile  = "Coeffs.txt"                   ### produced by Merope (coefficients)

fft.setNbOfThreads(number_of_processors)

def thermic_tmfft(output):
    g = fft.Grid(fft.VTKRead(voxellation_of_zones)) ### read the voxellation
    med = fft.Medium(g, "grdT", "Flx", "THERMAL")   ### experiment
    med.declareParamsT("lam")                       ### thermal conductivity

    # read the coefficient
    file_here = open(coeffFile,"r")
    lines = file_here.readlines()                   
    the_coeffs = [float(line) for line in lines]
    file_here.close()

    # set the coefficients
    for i, coeff in enumerate(the_coeffs):
        med[i].setConstant("lam",coeff)

    # Parametrization for thermal experiment
    slv = fft.TSolver(med)
    slv.setPrecision(1e-4)
    slv.setMaxIterations(500000)

    # Choose the reference medium for preconditioning
    lam1 = max(the_coeffs)
    lam0 = min(the_coeffs)
    slv.setLam0(sqrt(lam1*lam0))
    slv.Acceleration("ONEITER", 2)  # Acceleration factor of Milton and Eyre (2)
    #slv.setLam0(0.5*(lam0 + lam1))
    #slv.setAnderson("DF", 3)
    slv.setFiltering("ROTATION")
    
    # For vtk outputs
    VTKout = slv.setVTKout()
    VTKout.setBasename("thermique")
    VTKout.TField("T")
    VTKout.newField("Flx")
    VTKout.newField("lam")
    
    # Compute the equivalent conductivity
    slv.conductivityMatrix(293, 100, output)  # Reference temperature, gradient of temperature, output file
    file = open(output, 'a')
    file.write(str(g))
    file.close()

os.chdir(arch.resultFolder)
thermic_tmfft("tmfft_output.txt")
os.chdir("../")
