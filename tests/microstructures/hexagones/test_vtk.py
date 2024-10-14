# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
print(vtkreader_merope.compare_error("Zone.vtk","ref/Zone.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("Zone_bis.vtk","ref/Zone_bis.vtk", "MaterialId", 0.01))

