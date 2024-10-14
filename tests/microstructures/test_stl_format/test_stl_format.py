# -*- coding:utf8 -*-
#
# Copyright : see License.txt
#
# Tests stl_format

import sac_de_billes
import merope
import time

tic0 = time.time()

### Define of the micro-structure

sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength([10,10,10])
sphIncl.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[0.5,0.2],[0.25,0.2]], [1,2])

multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(sphIncl)

structure = merope.Structure_3D(multiInclusions)
gridParameters = merope.vox.create_grid_parameters_N_L_3D([100,100,100], [10,10,10])
my_printer = merope.vox.vtk_printer_3D()

coeffs = [i for i in range(max(multiInclusions.getAllIdentifiers()))]

### Center, phase
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)
grid.convert_to_stl_format()
internal_field = grid.get_as_list()
print("### Center, phase")
print(internal_field[0])

### Center, real
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)
grid.apply_coefficients(coeffs)
grid.convert_to_stl_format()
internal_field = grid.get_as_list()
print("### Center, real")
print(internal_field[0])

### Iso, phase
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Average)
grid.convert_to_stl_format()
internal_field = grid.get_as_list()
print("### Iso, phase")
print(internal_field[0])

### Iso, real
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Average)
grid.apply_coefficients(coeffs)
grid.convert_to_stl_format()
internal_field = grid.get_as_list()
print("### Iso, real")
print(internal_field[0])

### AnIso, phase
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Laminate)
grid.convert_to_stl_format()
internal_field = grid.get_as_list()
print("### AnIso, phase")
print(internal_field[0])

### AnIso, real
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Laminate)
grid.apply_coefficients(coeffs)
grid.convert_to_stl_format()
internal_field = grid.get_as_list()
print("### AnIso, real")
print(internal_field[0])