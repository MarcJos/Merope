# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 09/09/2022
#
# Copyright : see License.txt
#
# SpheroPolyhedron

import sac_de_billes
import merope
from math import *


L = [4., 4., 4.]
nbVox = [64, 64, 64]

#spheroPolyhedron = merope.microInclusion.SpheroPolyhedronFactory_3D.fromVertices(int(1), [[0., 0., 0.], [1., 0, 0.], [0, 1., 0], [0., 0., 1.]], [[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]], float(0.5))
A1 = [1.5, 1.5, 1.5]
A2 = [2.5, 1.5, 1.5]
A3 = [1.5, 2.5, 1.5]
A4 = [1.5, 1.5, 2.5]
faceVert =[[0, 2, 1], [0, 1, 3], [1, 2, 3], [0, 3, 2]]
minkowskiRadius = 0.5
spheroPolyhedron = merope.microInclusion.SpheroPolyhedronFactory_3D().fromVertices(int(1), [A1, A2, A3, A4], faceVert, minkowskiRadius)

multiInclusions = merope.MultiInclusions_3D()
multiInclusions.setInclusions([spheroPolyhedron], L)

grid = merope.Voxellation_3D(multiInclusions)
grid.setPureCoeffs([1,2])
grid.proceed(nbVox)
grid.printFile("Zone.vtk","Coeffs.txt")

