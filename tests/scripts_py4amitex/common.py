# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/09/2021
#
# Copyright : see License.txt
#
# Build a .geo file to be meshed by gmsh
# The geometry is made of polyhedra


import sac_de_billes
import merope
import numpy as np

nbSpheres = 16 
distMin = 0.05
randomSeed = 0

def create_structure(L):
    theSpheres = sac_de_billes.fillMaxRSA_3D(sac_de_billes.Tore, L, nbSpheres, randomSeed, distMin)
    sphInc = merope.LaguerreTess_3D(L,theSpheres)
    mi = merope.MultiInclusions_3D()
    mi.setInclusions(sphInc)
    mi.changePhase([i for i in mi.getAllIdentifiers() if i%2 == 0], [0 for i in mi.getAllIdentifiers() if i%2 == 0])
    mi.changePhase([i for i in mi.getAllIdentifiers() if i%2 == 1], [1 for i in mi.getAllIdentifiers() if i%2 == 1])
    structure = merope.Structure_3D(mi)
    return structure

def create_grid(structure, L, N3D, voxelRule):
    gridParameters = merope.vox.create_grid_parameters_N_L_3D(N3D, L)
    grid = merope.vox.GridRepresentation_3D(structure, gridParameters, voxelRule)
    return grid

def print_grid(structure, L, N3D):
    grid = create_grid(structure, L, N3D, merope.vox.VoxelRule.PolyGeom)
    grid.apply_coefficients(therm_pure_coeffs)
    grid.convert_to_Iso_format()
    grid.apply_homogRule(merope.HomogenizationRule.Voigt)
    my_printer = merope.vox.vtk_printer_3D()
    my_printer.printVTK(grid, "Zone.vtk")
    my_printer.printVTK_segmented(grid, "zoneIds.vtk", "Coeffs.txt")


def meca_pure_coeffs(iphase):
    if iphase == 0:
        return [100e6, 100e6]
    elif  iphase == 1:
        return [200e6, 200e6]
    else:
        raise Exceptiont("Unknown phase")

def meca_ref_coeffs(j):
    return 0.5 * (meca_pure_coeffs(0)[j] + meca_pure_coeffs(1)[j])

therm_pure_coeffs = [1, 2]
therm_ref_coeff = np.max(therm_pure_coeffs)

number_phases = 2
voxelRule = merope.vox.VoxelRule.Laminate
direction = [1, 0, 0]

if __name__=="__main__":
    L = [1, 1, 1]
    voxelRule = merope.vox.VoxelRule.Laminate
    GRID_DIMS = [4 for _ in range(3)]
    structure = create_structure(L)
    grid = create_grid(structure, L, GRID_DIMS, voxelRule)
    print_grid(structure, L, GRID_DIMS)
