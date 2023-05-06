#!/usr/bin/python
# -*- coding: utf-8 -*- 
# 
# benchmark.py
# Marc Josien
# 01/10/2021
#
# Compare performances of tmfft and mérope
#
# Copyright : see License.txt

import benchmark as benchmark

list_time = []

L0 = 20
nVox0 = 32
nbThreads = 4

jmin = 0
imin = 0
jmax = 5
imax = 6

with open("Time.txt", "w") as fic:
    fic.write("nVox0 L0 nbSpheres nbThreads timeMerope timeTMFFT \n")

for j in range(jmin,jmax):
    nVox = 2**j * nVox0
    for i in range(imin,imax):
        L = 2**i * L0
        result = benchmark.testVoxel(nbThreads, L, nVox, False)
        times = result.timeList
        nbSpheres = result.nbSpheres
        list_time.append(times)
        with open("Time.txt", "a") as fic:
            fic.write(str(nVox) + " ")
            fic.write(str(L) + " ")
            fic.write(str(nbSpheres) + " ")
            fic.write(str(nbThreads) + " ")
            fic.write(str(times[0]) + " " + str(times[1]) + "\n")

