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


import sys
import os
sys.path.append(os.path.dirname(os.path.abspath(__file__)))

import JDD_1 as jdd1
import JDD_2 as jdd2
import JDD_3 as jdd3
import JDD_4 as jdd4




nbVox = [100, 100, 100]

def print_spheroPoly(L, vertices, faces_vertices_indexes, minkowskiRadius, nameFileVTK, nbVox):
    spheroPolyhedronFactory = merope.microInclusion.SpheroPolyhedronFactory_3D()
    spheroPolyhedron = spheroPolyhedronFactory.fromVertices(1, vertices, faces_vertices_indexes, minkowskiRadius)
    spheroPolyhedronList = [spheroPolyhedron]

    spheroPolyInclusions = merope.SpheroPolyInclusions_3D()
    spheroPolyInclusions.setLength(L)
    spheroPolyInclusions.setInclusions(spheroPolyhedronList)
    
    multiInclusions = merope.MultiInclusions_3D()
    multiInclusions.setInclusions(spheroPolyInclusions)
    multiInclusions.setMatrixPhase(0)

    structure = merope.Structure_3D(multiInclusions)

    voxellation = merope.Voxellation_3D(structure)
    voxellation.proceed(nbVox)
    voxellation.printFile(nameFileVTK,"Coeffs.txt")

print_spheroPoly(jdd1.L, jdd1.vertices, jdd1.faces_vertices_indexes, jdd1.minkowskiRadius, "jdd1.vtk", nbVox)
print_spheroPoly(jdd2.L, jdd2.vertices, jdd2.faces_vertices_indexes, jdd2.minkowskiRadius, "jdd2.vtk", nbVox)
print_spheroPoly(jdd3.L, jdd3.vertices, jdd3.faces_vertices_indexes, jdd3.minkowskiRadius, "jdd3.vtk", nbVox)
print_spheroPoly(jdd4.L, jdd4.vertices, jdd4.faces_vertices_indexes, jdd4.minkowskiRadius, "jdd4.vtk", nbVox)
