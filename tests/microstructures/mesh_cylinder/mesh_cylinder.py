# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 03/10/2024
# Copyright : see License.txt
#
# Build a .geo file to be meshed by gmsh
# The geometry is made of cylinders

import sac_de_billes
import merope
import time
import random
import numpy as np

L = [5, 5, 5]

seed = 0
mindist = 0.2

theSpheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seed, [[1, 1]], [1], mindist)

def create_cylinder_inc(pt1, pt2, radius):
    axis = merope.geometry.Segment_3D([pt1, pt2])
    cylinder = merope.geometry.Cylinder(axis, radius)
    cylinderInc = merope.microInclusion.CylinderInclusions(cylinder)
    return cylinderInc

radius = 0.4
height = np.sqrt(1-radius*radius)

random.seed(seed)
list_cylinders = []

for s in theSpheres:
    orient = [random.random() for i in range(3)]
    norm = np.sqrt(orient[0] * orient[0] + orient[1] * orient[1] + orient[2] * orient[2])
    orient = [orient[i] / norm for i in range(3)]
    pt1 = [s.center[i] + 0.5 * height * orient[i] for i in range(3)]
    pt2 = [s.center[i] - 0.5 * height * orient[i] for i in range(3)]
    list_cylinders.append(create_cylinder_inc(pt1, pt2, radius))

cylinderInclusions = merope.CylinderInclusions_3D()
cylinderInclusions.setLength(L)
cylinderInclusions.setInclusions(list_cylinders)

multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions(cylinderInclusions)
multiInclusions.changePhase(multiInclusions.getAllIdentifiers(), [2 for i in multiInclusions.getAllIdentifiers()])
multiInclusions.setMatrixPhase(1)
multiInclusions.addLayer(multiInclusions.getAllIdentifiers(), 3, 0.15)
multiInclusions.addLayer(multiInclusions.getAllIdentifiers(), 4, 0.15)

### voxellize it

structure = merope.Structure_3D(multiInclusions)
gridParameters = merope.vox.create_grid_parameters_N_L_3D([256,256,256], L)
grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(grid, "Zone.vtk")

### mesh it

meshGenerator = merope.mesh.MeshGenerator()
meshGenerator.setMeshOrder(2)
meshGenerator.setMeshSize(0.5 * mindist)
meshGenerator.setMultiInclusions(multiInclusions)
meshGenerator.do_not_mesh([2])
meshGenerator.set_nameOutput(["cylinders.vtk"])
meshGenerator.write("mesh_cylinders.geo")
