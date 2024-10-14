# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope

liste_names = ["orig_crystal.vtk", "smaller_crystal.vtk", "larger_smaller_crystal_0.vtk", "larger_smaller_crystal.vtk", "orig_sphere.vtk", "smaller_sphere.vtk", "larger_smaller_sphere_0.vtk", "larger_smaller_sphere.vtk"]

for name in liste_names:
    print(vtkreader_merope.compare_error(name,"ref/" + name , "MaterialId", 0.01))

