# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/07/2024
#
# Copyright : see License.txt
#
# Create a configuration of penny-shape inclusions
# All are inside non-intersecting spheres
# All are oriented in the same direction

import sac_de_billes
import merope
import time
import numpy as np


tic = time.time()

l_0 = 1

ratio = 0.01 # height/diameter
L = [l_0, l_0,  ratio *l_0]
cylinder_radius = 0.02
cylinder_height = 2 * ratio * cylinder_radius
seed = 0

# In the dilated world, cylinders are as large as high
L_dilated_by_ratio = [L[0], L[1], L[2]/ratio]
new_radius = np.sqrt(2) * cylinder_radius
theCenters_dilated = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L_dilated_by_ratio, seed, [[new_radius, 1]], [0], 0)

# back to normal world
for s in theCenters_dilated:
    s.center = [s.center[0], s.center[1], s.center[2] * ratio]

theCenters = theCenters_dilated

def fromSphere_to_cylinder(sphere):
    center = sphere.center
    axis = merope.geometry.Segment_3D([[center[0], center[1], center[2] - 0.5 * cylinder_height], [center[0], center[1], center[2] + 0.5 * cylinder_height]])
    cylinder = merope.geometry.Cylinder(axis, cylinder_radius)
    return cylinder


theCylinders = [fromSphere_to_cylinder(s) for s in theCenters]

cylinderInclusions = merope.CylinderInclusions_3D()
cylinderInclusions.setLength(L)
cylinderInclusions.setInclusions([merope.microInclusion.CylinderInclusions(cylinder) for cylinder in theCylinders])

multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(cylinderInclusions)
multiInclusions.setMatrixPhase(1)


structure = merope.Structure_3D(multiInclusions)

gridParameters = merope.vox.create_grid_parameters_N_L_3D([300,300,300], L)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(grid, "Zone.vtk")


print(time.time()- tic)