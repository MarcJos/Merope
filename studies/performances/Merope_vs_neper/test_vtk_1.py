# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021
# Copyright : see License.txt

import vtkreader_merope
percentageDifferentVoxels = vtkreader_merope.compareVoxellationPercentage("Merope128_128_10000.vtk","Neper128_128_10000.vtk", "MaterialId")
print(percentageDifferentVoxels)
#print(vtkreader_merope.compare_error("Merope128_128_10000.vtk","Neper128_128_10000.vtk", "MaterialId", 1))

