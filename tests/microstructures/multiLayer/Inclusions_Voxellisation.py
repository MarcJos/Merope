# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 10/06/2021
#
# Copyright : see License.txt

import merope
import sac_de_billes
import time

tic0 = time.time()

n3D = 256

L = [10,10,10]

sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength(L)
sphIncl.fromHisto(0, sac_de_billes.TypeAlgo.RSA, 0., [[0.5,0.2],[0.25,0.2]], [1,2])

mIncl = merope.MultiInclusions_3D()
mIncl.setInclusions(sphIncl)


list0 = mIncl.getAllIdentifiers()
list1 = mIncl.getIdentifiers([1,2])
list2 = mIncl.getIdentifiers([1])

mIncl.addLayer(list1,3, 0.1)
mIncl.addLayer(list2,4, 0.1)

list3 = mIncl.getIdentifiers([1])


###############################################
### Test coherence
list1.sort()
list0.sort()
if len(list0) != len(list1):
    raise NameError("Problem0")
for i, n in enumerate(list0):
    if list0[i] != list1[i]:
        raise NameError("Problem1")
###


### Test coherence entre 2 instants différents
list2.sort()
list3.sort()
if len(list2) != len(list3):
    raise NameError("Problem2")
for i, n in enumerate(list2):
    if list2[i] != list3[i]:
        raise NameError("Problem3")
###
###############################################

structure = merope.Structure_3D(mIncl)

gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D, n3D, n3D], L)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(grid, "Zone_Incl.vtk")

print(time.time() - tic0)
