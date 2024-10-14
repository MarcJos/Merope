# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
print(vtkreader_merope.compare_error("GaussianField.vtk","ref/GaussianField.vtk", "MaterialId", 4))

