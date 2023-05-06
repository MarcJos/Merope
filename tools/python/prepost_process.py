#!/usr/bin/python
# coding: utf-8
# 2021-10
# Prepocressing and post-processing functions
#
# Copyright : see License.txt

import csv
import sac_de_billes


def readSpheres(file):
    theSpheres = []
    with open(file,"r") as inputFic:
        reader = csv.reader(inputFic, delimiter=' ')
        for row in reader:
            center = [float(row[i]) for i in range(0,3)]
            radius = float(row[3])
            if len(row) == 4:
                phase = 0
            else :
                phase = int(row[4])
            theSpheres.append(sac_de_billes.Sphere_3D(center, radius, phase))
    return theSpheres


def readDiffusionCoeffs(coeffFileName):
    result = []
    with open(coeffFileName, "r") as fic:
        for line in fic:
            result.append(float(line))
    return result
