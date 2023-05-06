# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
refs = ["jdd1_ref.vtk", "jdd2_ref.vtk", "jdd3_ref.vtk", "jdd4_ref.vtk"]
simu = ["jdd1.vtk", "jdd2.vtk", "jdd3.vtk", "jdd4.vtk"]
for i in range(len(refs)):
    print(vtkreader_merope.compare_error(refs[i], simu[i], "MaterialId", 1))

