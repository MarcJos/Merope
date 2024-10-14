# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
print(vtkreader_merope.compare_error("Zone_struct_1.vtk","ref/Zone_struct_1.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("Zone_struct_2.vtk","ref/Zone_struct_2.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("Zone_struct_mask.vtk","ref/Zone_struct_mask.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("Zone_struct_finale.vtk","ref/Zone_struct_finale.vtk", "MaterialId", 0.01))

