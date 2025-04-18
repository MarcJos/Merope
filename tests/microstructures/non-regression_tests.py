# -*- coding:utf8 -*-
#
# Copyright : see License.txt
#
# Non-regression tests

import os
import time
import sys
import importlib
import subprocess

def go_to_dir(name_dir):
    time.sleep(1)
    os.chdir(name_dir)
    print(name_dir)
    time.sleep(1)

def execute_tests(name_dir):
    go_to_dir(name_dir)
    subprocess.run(["python3", name_dir + ".py"], check=True)
    importlib.import_module(name_dir + ".test_vtk")
    go_to_dir("../")

list_of_test_names = [
    'buildVoxellation',
    'closestNeighbors',
    'coated_inclusions', 
    'Cpp_tests',
    'cylinder',
    'enlarge',
    'fields_output',
    'hexagones', 
    'inclusions', 
    'intersectingSpheres',
    'laguerre_tess_non_periodic',
    'largeSphere',
    'many_cylinders', 
    'many_cylinders_random_orient',
    'mesh_cylinder',
    'mesh_poly_0',
    'mesh_poly_1',
    'mesh_spheres_0',
    'mesh_spheres_1',
    'mesh_spheres_2',
    'multiLayer', 
    'optimize_Laguerre_2D', 
    'optimize_Laguerre_tess',
    'mask_planeGeom',
    'planeGeom_check',
    'polyCrystal', 
    'polyCrystal_2D', 
    'polyCrystal_filamentaire',
    'recursive_structure', 
    'spheroPolyhedron', 
    'test_stl_format',
    'variousSpheroPolyhedron',
    'voro_wall']

for test_name in list_of_test_names:
    execute_tests(test_name)



