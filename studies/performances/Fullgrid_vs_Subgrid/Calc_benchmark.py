#!/usr/bin/python
# -*- coding: utf-8 -*- 
# Author: M. Josien
#
# Copyright : see License.txt

from benchmark import *

#os.system("export OMP_NUM_THREADS=4")
#Calculs sur 4 threads pour Neper

TimeList_3D = []
TimeList_3D_Slice = []
List_L0 = [32,64,128,256,512,1024]
List_nSpheres = [10,100,1000,10000]
FichierSortie = './TimeVox.res'
off = open(FichierSortie,'w')
off.write("nbVox0  nSpheres Time_3D_full Time_3D_SubGrid\n")


for L0 in List_L0:
    nVox0 = L0
    for nSpheres in List_nSpheres:
        Test = testVoxel(L0, nVox0, nSpheres)
        TimeList_3D.append(Test.timeList[0])
        TimeList_3D_Slice.append(Test.timeList[1])
        off.write("\n"+"%.13e"% L0+" "+"%.13e"% nSpheres+" "+"%.13e"% Test.timeList[0]+" "+"%.13e"% Test.timeList[1])
    off.write("\n")


off.close()

