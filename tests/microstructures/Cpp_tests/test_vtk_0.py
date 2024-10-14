# -*- coding:utf8 -*-
#
# Test the vtk
# Author: M. Josien
# Date: 07/2021

import vtkreader_merope
liste_name = ["inclusions.vtk", "poly3D-gdVolume.vtk", "poro2D_1.vtk", "poro2D_2.vtk", "poro3D_2_0.vtk", "poro3D_2.vtk", "poro3D-gdVolume_2.vtk", "poro3D_spheres_2_sym.vtk", "poro3D_spheres_2.vtk", "Spheres.vtk", "field_struc_gauss.vtk", "struc_incl.vtk", "field_struc_mask.vtk", "totalStruct_real.vtk", "Zone_Crystal_extraction.vtk", "Zone_extraction.vtk", "Zone_Inclusions_extraction.vtk", "poro3D_1__.vtk", "textured.vtk"]

liste_name_ref = ["ref_" + name for name in liste_name]

import os

for i in range(len(liste_name)):
    print(liste_name[i])
    print(liste_name_ref[i])
    print(vtkreader_merope.compare_error_double(liste_name[i], liste_name_ref[i], "MaterialId", 0.01))



