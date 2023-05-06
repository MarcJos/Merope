#!/usr/bin/python
# -*- coding: utf-8 -*- 
#
# Author: M. Josien
#
# Copyright : see License.txt
#

from benchmark import *
from os import system
os.system("export OMP_NUM_THREADS=1")
#Calculs sur 1 thread pour Neper

os.system("echo $OMP_NUM_THREADS")

TimeList_Neper = []
TimeList_Merope = []
List_L0 = [32,64,128,256, 512]
List_nSpheres = [10,100,1000,10000, 100000]
FichierSortie = './TimeVox.res'
off = open(FichierSortie,'w')
off.write("#nSpheres L0(=nVox0) Time_Neper Time_Merope\n")


for nSpheres in List_nSpheres:
    for L0 in List_L0:
        nVox0 = L0
        Test = testVoxel(L0, nVox0, nSpheres)
        TimeList_Neper.append(Test.timeList[0])
        TimeList_Merope.append(Test.timeList[1])
        off.write("\n"+"%.13e"% nSpheres+" "+"%.13e"% L0+" "+"%.13e"% Test.timeList[0]+" "+"%.13e"% Test.timeList[1])
#    off.write("\n")


off.close()

