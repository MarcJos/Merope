# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 29/08/2022
#
# Copyright : see License.txt
#
# Get a phase field from a microstructure

# Import Merope & tmfft
import sac_de_billes
import merope
import archi_merope as arch

### Define of the micro-structure

L = [10, 10, 10]
sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength(L)
sphIncl.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[0.5,0.2],[0.25,0.2]], [1,2])

### Print outputs

sphIncl.printDump("inclusions.dump")
sphIncl.printPos("inclusions.pos")
sphIncl.printCSV("inclusions.csv")
sphIncl.printFracVol("volume_fractions.txt")

### Build the voxellation

multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(sphIncl)

nbVox = [100,100,100]

### Play with the phaseField
def printHashtags():
  for i in range(10):
    print("##############")


structure = merope.Structure_3D(multiInclusions)

gridParameters = merope.vox.create_grid_parameters_N_L_3D(nbVox, L)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Average)

copy_grid = merope.vox.GridRepresentation_3D(grid) #deep copy
copy_grid.convert_to_stl_format()
isoPhaseField = copy_grid.get_as_list()

printHashtags()
expected_list = [[(2, 1.0)], [(2, 1.0)], [(2, 1.0)], [(2, 0.5141177757852943), (0, 0.4858822242147057)], [(2, 0.4045971926096369), (0, 0.5954028073903631)]]
ii = 0
for i in range(3,8):
  j = 5
  k = 12
  print("The content of the voxel at position in discrete coordinates [i, j ,k] = " + str([i,j,k]) + " is ")
  linear_discrete_coordinate = merope.vox.get_linear_index_3D([i, j, k], nbVox)
  current_phase = isoPhaseField[linear_discrete_coordinate]
  print(current_phase)
  for l in range(len(expected_list[ii])):
    if current_phase[l][0] != expected_list[ii][l][0] or abs(current_phase[l][0] - expected_list[ii][l][0]) > 1e-6:
      raise Exception("Unexpected")
  print("which is a list of [phase, volumeFraction]")
  ii += 1
printHashtags()


grid.apply_homogRule(merope.HomogenizationRule.Voigt, [1,2,3])

numpy_converter = merope.vox.NumpyConverter_3D()
realField = numpy_converter.compute_RealField(grid)



printHashtags()
for i in range(3,8):
  j = 5
  k = 12
  print("The content of the voxel at position in discrete coordinates [i, j ,k] = " + str([i,j,k]) + " is ")
  linear_discrete_coordinate = merope.vox.get_linear_index_3D([i, j, k], nbVox)
  print(realField[linear_discrete_coordinate])
  print("which is a list of [phase, volumeFraction]")
printHashtags()

### Print the grid
my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK_segmented(grid, "Zone.vtk","Coeffs.txt")

