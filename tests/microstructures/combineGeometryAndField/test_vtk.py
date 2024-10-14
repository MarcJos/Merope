# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
print(vtkreader_merope.compare_error("CombineField.vtk","ref/CombineField.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("Field0.vtk","ref/Field0.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("Field1.vtk","ref/Field1.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("Field2.vtk","ref/Field2.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("Inclusions.vtk","ref/Inclusions.vtk", "MaterialId", 0.01))

