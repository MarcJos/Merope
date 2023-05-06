#!/usr/bin/python
# coding: utf-8
# 07/2020
#
# Computes the homogenized matrix with outputs of Mérope
# Mainly for comparing with amitex
#
# Copyright : see License.txt

import tmfft as fft 
from math import *
import numpy as np

magic_constant_nabla_T = 100
default_temperature = 293

def tmfft_thermo(vtkFile, the_coeffs, NB_OF_PROC, output, see_output = False,  T_0 = default_temperature, precision=1e-4):
    ### How many processors?
    fft.setNbOfThreads(NB_OF_PROC)
    
    g = fft.Grid(fft.VTKRead(vtkFile))
    # Création du milieu pour la thermique
    med = fft.Medium(g, "grdT", "Flx", "THERMAL")
    # Conductivité thermique
    med.declareParamsT("lam")

    for i, coeff in enumerate(the_coeffs):
        med[i].setConstant("lam",coeff)

    # Solveur pour la thermique
    slv = fft.TSolver(med)
    slv.setPrecision(precision)
    slv.setMaxIterations(500000)

    # Paramètres des accélérations
    lam1 = max(the_coeffs)
    lam0 = min(the_coeffs)
    slv.setLam0(sqrt(lam1*lam0))
    slv.Acceleration("ONEITER", 2)  # Facteur d'acceleration Milton et Eyre (2)
#    slv.setFiltering("ROTATION")
    
    # Pour les sorties VTK
    if see_output:
        VTKout = slv.setVTKout()
        VTKout.setBasename("thermique")
        VTKout.TField("T")
        VTKout.newField("Flx")
        VTKout.newField("lam")
    print(g)
    
    # Conductivité équivalente
    slv.conductivityMatrix(T_0, magic_constant_nabla_T, output)  # Temperature de référence, gradient de temperature, nom du fichier de sortie
    file = open(output, 'a')
    file.write(str(g))
    file.close()
