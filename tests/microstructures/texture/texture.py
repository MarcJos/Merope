# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 20/09/2024
#
# Copyright : see License.txt
#
# Build textured microstructure 

import sac_de_billes
import merope
import random
import numpy as np
from lambda_function_builder import *
from numba import cfunc


seed = 0  # Random number generations seed for seperating spheres positionning
l3D = 1 # RVE dimensions
L = [ l3D, l3D, l3D ]
desiredRPhi = [[0.04,  0.1]]
theSpheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.RSA, sac_de_billes.NameShape.Tore, L, seed, desiredRPhi, [1], 0.)

lag = merope.LaguerreTess_3D(L, theSpheres)
lag.setAspRatio([2, 1., 1.])
mi = merope.MultiInclusions_3D()
mi.setInclusions(lag)
identifiers = mi.getAllIdentifiers()
mi.changePhase(identifiers, identifiers)

structure = merope.Structure_3D(mi)


# texturing
# building the map
def pickOnSphere_3D():
    while True:
        x = [random.random() for i in range(3)]
        norme = np.sqrt(x[0]*x[0] + x[1]*x[1] + x[2]*x[2])
        if norme < 1:
            return [x[i]/norme for i in range(3)]

random.seed(seed)
identifier_to_orientation = { id : pickOnSphere_3D() for id in identifiers}
# texturer
period = 0.005

# Define a large array in the main code
# Because dictionnary are not supported inside a cfunc
max_index = max(identifiers) + 1
d = np.zeros((max_index,3), dtype=np.float64)
for k, v in identifier_to_orientation.items():
    d[k][0] = v[0]
    d[k][1] = v[1]
    d[k][2] = v[2]

#@cfunc(types.float64(types.CPointer(types.float64), types.int64))
def texturer(x, phase):
    z = d[phase]
    return np.sin((x[0] * z[0]  + x[1] * z[1]  + x[2] * z[2]) / period)
merope_texturer = texture_py2c(texturer, 3)

#! voxellize basic structure
n = 128
nbNodes = [n, n, n]

grid_parameters = merope.vox.create_grid_parameters_N_L_3D(nbNodes, L)
phaseGrid = merope.vox.GridRepresentation_3D(structure, grid_parameters, merope.vox.VoxelRule.Center)

phaseGrid.apply_texture(merope_texturer.get_funcPointer())

my_printer = merope.vox.vtk_printer_3D()
my_printer.printVTK(phaseGrid, "phases.vtk")