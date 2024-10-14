# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
print(vtkreader_merope.compare_error("result_folder/Phases.vtk","ref/Phases.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("result_folder/Composite.vtk","ref/Composite.vtk", "MaterialId", 0.01))

