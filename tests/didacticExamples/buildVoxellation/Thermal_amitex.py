# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 27/09/2021
#
# Copyright : see License.txt
#
# Launches AMITEX computation
# 
import os
import archi_merope as arch
import interface_amitex_fftp.amitex_wrapper as amitex
import interface_amitex_fftp.post_processing as amitex_out


os.chdir(arch.resultFolder)

### Launch amitex
number_of_processors = 2                   ### for parallel computing
voxellation_of_zones = "Composite.vtk"      ### produced by Merope
### Implicitly uses the file "Coeffs.txt" to compute the thermal matrix
amitex.computeThermalCoeff(voxellation_of_zones, number_of_processors)
### See the homogenized matrix
homogenized_matrix = amitex_out.printThermalCoeff(".")


os.chdir("../")
