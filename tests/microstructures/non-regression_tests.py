# -*- coding:utf8 -*-
#
# Copyright : see License.txt
#
# Non-regression tests

import os
import time
import sys

def go_to_dir(name_dir):
#    sys.path.insert(0, os.getcwd() + "/" + name_dir)
#    print(os.getcwd() +  name_dir)
    time.sleep(1)
    os.chdir(name_dir)
    print(name_dir)
    time.sleep(1)

go_to_dir("variousSpheroPolyhedron")
import variousSpheroPolyhedron.variousSpheroPolyhedron
import variousSpheroPolyhedron.test_vtk
go_to_dir("../")

go_to_dir("recursive_structure")
import recursive_structure.recursive_structure
go_to_dir("../")


go_to_dir("polyCrystal_filamentaire")
import polyCrystal_filamentaire.polyCrystal_filamentaire
go_to_dir("../")

go_to_dir("polyCrystal_2D")
import polyCrystal_2D.polyCrystal_2D
go_to_dir("../")

go_to_dir("polyCrystal")
import polyCrystal.polyCrystal
go_to_dir("../")

go_to_dir("multiLayer")
import multiLayer.polyCrystal_Voigt
import multiLayer.Inclusions_Voxellisation
go_to_dir("../")

go_to_dir("inclusions")
import inclusions.Inclusions_Voxellisation
go_to_dir("../")

go_to_dir("hexagones")
import hexagones.hexagones
go_to_dir("../")

go_to_dir("coated_inclusions")
import coated_inclusions.coated_inclusions
go_to_dir("../")

go_to_dir("gaussian")
import gaussian.gaussian
import gaussian.test_vtk
go_to_dir("../")

go_to_dir("prescribedField")
import prescribedField.prescribedField
go_to_dir("../")


go_to_dir("Cpp_tests")
import Cpp_tests.Cpp_tests
go_to_dir("../")

go_to_dir("combineGeometryAndField")
import combineGeometryAndField.combineGeometryAndField
go_to_dir("../")


go_to_dir("optimize_Laguerre_tess")
import optimize_Laguerre_tess.optimize_Laguerre_tess
go_to_dir("../")

