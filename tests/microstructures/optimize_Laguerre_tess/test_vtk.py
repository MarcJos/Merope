# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
print(vtkreader_merope.compare_error_double("all_tessels_with_same_volume.vtk","ref/all_tessels_with_same_volume.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error_double("tessels_with_weight.vtk","ref/tessels_with_weight.vtk", "MaterialId", 0.15))
print(vtkreader_merope.compare_error_double("original_tessellation.vtk","ref/original_tessellation.vtk", "MaterialId", 0.15))
