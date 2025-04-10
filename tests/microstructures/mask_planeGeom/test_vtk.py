# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope

nameFileList = ["average.vtk", "core.vtk", "layer.vtk", "mask.vtk", "planeGeom.vtk"]

for nameFile in nameFileList:
    print("Compare : " + nameFile)
    print(vtkreader_merope.compare_error(nameFile,"ref/" + nameFile, "MaterialId", 0.01))



