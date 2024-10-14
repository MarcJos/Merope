# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
print(vtkreader_merope.compare_error("Zone_Voigt.vtk","ref/Zone_Voigt.vtk", "MaterialId", 0.01))
print(vtkreader_merope.compare_error("Zone_Incl.vtk","ref/Zone_Incl.vtk", "MaterialId", 0.01))

