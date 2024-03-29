import sac_de_billes
import merope
from math import *
import ctypes
from numba import cfunc, types
import numpy as np
import pybind11 as pyb

import time

def pb(nb_proc):
    my_time = time.time()
    merope.setNbOfThreads(nb_proc) 
    
    
    n3D = 64
    L = [10, 10, 10]

    def printField(field, L, n3D, name_vtk):
      cartesianField = merope.CartesianField_3D(field, L)
      structureG = merope.FieldStructure_3D(cartesianField)
      vox = merope.Voxellation_3D(structureG)
      vox.proceed([n3D, n3D, n3D])  
      vox.printFieldFile(name_vtk)
    
    
    def createGaussianField(lambda_multiplication, name, seed):
      covariance = lambda x : exp((-x[0]**2 - x[1]**2 - x[2]**2) * lambda_multiplication)
      signature_cov = types.float64(types.CPointer(types.float64)) # double(double*)
      c_cov = cfunc(signature_cov, nopython=True)(covariance)
    
      nonLin = lambda x : lambda_multiplication *  exp(-x**2)
      signature_nl = types.float64(types.float64) # double(double)
      c_nln = cfunc(signature_nl, nopython=True)(nonLin)

      pt_adress_cov = merope.gaussianField.Interf_FuncPointer(c_cov.address, True)
      pt_adress_nln = merope.gaussianField.Interf_FuncPointer(c_nln.address, True)
    
      gaussianne = merope.gaussianField.GaussianField_3D(pt_adress_cov, pt_adress_nln)
      gaussianne.seed = seed
      
      ### print the desired covariance
      covField = merope.ScalarField_3D(covariance)
      printField(covField, L, n3D, name + "desired_cov.vtk")
      ### print the actual covariance
      #  seenCovField =  merope.gaussianField.NumericalCovariance_3D(covariance)
      #  printField(seenCovField, L, n3D, name + "real_cov.vtk")
    
      cartesianField = merope.CartesianField_3D(gaussianne, L)
      structureG = merope.FieldStructure_3D(cartesianField)
      vox = merope.Voxellation_3D(structureG)
      vox.proceed([n3D, n3D, n3D])  
      vox.printFile(name + "GaussianField.vtk", name + "Coeffs.txt")  
      ### get the Gaussian Field as an numpy.array
      gaussianField = vox.getField()
      return structureG
      
    structureGaussianField_0 = createGaussianField(1, "", 1)
    structureGaussianField_1 = createGaussianField(5, "1_", 0)
    

    print("##################")
    print(time.time() - my_time)
    print("##################")



pb(1)

pb(2)
    
pb(4)    
    
    
    
    
    
