import merope
from math import *
from lambda_function_builder import *


import time

def pb(nb_proc):
    my_time = time.time()
    merope.setNbOfThreads(nb_proc) 
    
    
    n3D = 64
    L = [10, 10, 10]

    def printField(field, L, n3D, name_vtk):
      cartesianField = merope.vox.CartesianField_3D(field, L)
      structureG = merope.FieldStructure_3D(cartesianField)

      gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D, n3D, n3D], L)
      grid = merope.vox.GridRepresentation_3D(structureG, gridParameters, merope.vox.VoxelRule.Center)
      
      my_printer = merope.vox.vtk_printer_3D()
      my_printer.printVTK(grid, name_vtk)
    
    
    def createGaussianField(lambda_multiplication, name, seed):
      covariance = lambda x : exp((-x[0]**2 - x[1]**2 - x[2]**2) * lambda_multiplication)
      covariance_pointer = function_py2c(covariance, dim_from=3)
      interf_covariance = covariance_pointer.get_funcPointer()

      nonLin = lambda x : lambda_multiplication *  exp(-x**2)
      nonLin_pointer = function_py2c(nonLin, dim_from=1)
      interf_nln = nonLin_pointer.get_funcPointer()
      

      gaussianne = merope.gaussianField.GaussianField_3D(interf_covariance, interf_nln)
      gaussianne.seed = seed
      
      ### print the desired covariance
      covField = merope.ScalarField_3D(interf_covariance)
      printField(covField, L, n3D, name + "desired_cov.vtk")
      ### print the actual covariance
      seenCovField =  merope.gaussianField.NumericalCovariance_3D(interf_covariance)
      printField(seenCovField, L, n3D, name + "real_cov.vtk")
    
      cartesianField = merope.vox.CartesianField_3D(gaussianne, L)
      structureG = merope.FieldStructure_3D(cartesianField)

      gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D, n3D, n3D], L)
      grid = merope.vox.GridRepresentation_3D(structureG, gridParameters, merope.vox.VoxelRule.Center)
      
      my_printer = merope.vox.vtk_printer_3D()
      my_printer.printVTK_segmented(grid, name + "GaussianField.vtk", name + "Coeffs.txt")  
      ### get the Gaussian Field as a Mérope CartesianGrid
      gaussianField = grid.get_PureRealField()
      ### get the Gaussian Field as a numpy array
      numpy_converter = merope.vox.NumpyConverter_3D()
      numpy_filed = numpy_converter.compute_RealField(grid)
      return structureG
      
    structureGaussianField_0 = createGaussianField(1, "", 1)
    structureGaussianField_1 = createGaussianField(5, "1_", 0)
    

    print("##################")
    print(time.time() - my_time)
    print("##################")



pb(1)

pb(2)
    
pb(4)    

pb(8)    
    
    
    
    
    
