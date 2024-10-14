# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
print(vtkreader_merope.compare_error("Zone_struct_finale.vtk","ref/Zone_struct_finale.vtk", "MaterialId", 0.01))

