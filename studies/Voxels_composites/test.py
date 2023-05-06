# -*- coding:utf8 -*-
#
# Example for the use of composite voxels.
# We want to compute the conductivity of lead spheres coated with gold inside water.
#
# Author: M. Josien
# Date:  2021-12-16
#
#
# Copyright : see License.txt

from buildVoxellation import *
from parametrization import *
import merope
import os


multiInclusions = geometry() 

NList_subset = NList

for n in NList_subset:
    for i in range(0, len(voxelRuleList)):
        voxRule = voxelRuleList[i]
        voxName = voxelRuleNames[i]
        for j in range(0, len(homogRuleList)):
            hgRule = homogRuleList[j]
            hgName = homogRuleNames[j]
            if(i == 1 or j == 0): ## for center, only 1 name
                fileName = "Result" + voxName + "_" + hgName + ".res"
                coeff = wholeProcedure(n, voxRule, hgRule, multiInclusions)
                with open(fileName, "a") as fic:
                    fic.write(str(n) + " ")
                    fic.write(str(coeff))
                    fic.write("\n")
