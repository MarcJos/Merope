# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 24/09/2024
#
# Copyright : see License.txt
#
# Enlarge inclusions

import sac_de_billes
import merope

L = [1, 1, 1]
n3D =64
gridParameters = merope.vox.create_grid_parameters_N_L_3D([n3D,n3D,n3D], L)

seed = 0
spheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seed, [[0.2, 1]], [1], 0)

def print_multi(multi, nameFile):
    structure = merope.Structure_3D(multi)
    grid = merope.vox.GridRepresentation_3D(structure, gridParameters, merope.vox.VoxelRule.Center)
    my_printer = merope.vox.vtk_printer_3D()
    my_printer.printVTK(grid, nameFile)

###

polyCrystal = merope.LaguerreTess_3D(L, spheres)

sphIncl = merope.SphereInclusions_3D()
sphIncl.setLength(L)
sphIncl.setSpheres(spheres)

def make_XP(inclusions, name):
    multiInclusions = merope.MultiInclusions_3D()
    multiInclusions.setInclusions(inclusions)
    multiInclusions.changePhase(multiInclusions.getAllIdentifiers(), [1 for i in multiInclusions.getAllIdentifiers()])
    multiInclusions.addLayer(multiInclusions.getAllIdentifiers(), 2, 0.1)
    print_multi(multiInclusions, "orig_" + name + ".vtk")
    multiInclusions.enlarge(multiInclusions.getAllIdentifiers(), -0.05)
    print_multi(multiInclusions, "smaller_" + name + ".vtk")
    multiInclusions.enlarge(multiInclusions.getAllIdentifiers(), 0.025)
    print_multi(multiInclusions, "larger_smaller_" + name + "_0.vtk")
    multiInclusions.enlarge(multiInclusions.getAllIdentifiers(), 0.025)
    print_multi(multiInclusions, "larger_smaller_" + name + ".vtk")

make_XP(polyCrystal.toPolyInclusions(), "crystal")
make_XP(sphIncl, "sphere")
