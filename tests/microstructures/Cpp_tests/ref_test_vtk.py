# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
liste_name = ["gauss.vtk", "inclusions.vtk", "mask.vtk", "poly3D-gdVolume.vtk", "poro2D_1.vtk", "poro2D_2.vtk", "poro3D_1.vtk", "poro3D_2_0.vtk", "poro3D_2.vtk", "poro3D-gdVolume.vtk", "poro3D_spheres_2_sym.vtk", "poro3D_spheres_2.vtk", "Spheres.vtk", "struc_gauss.vtk", "struc_incl.vtk", "struc_mask.vtk", "test_vtk.py", "totalStruct.vtk", "Zone_Crystal_extraction.vtk", "Zone_extraction.vtk", "Zone_Inclusions_extraction.vtk"]

liste_name_ref = ["ref_" + name for name in liste_name]

import os

for i in range(len(liste_name)):
    os.system("cp " + liste_name[i] + " " + liste_name_ref[i])



