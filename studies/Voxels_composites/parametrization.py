# -*- coding:utf8 -*-
#
# Example for the use of composite voxels.
# We want to compute the conductivity of lead spheres coated with gold inside water.
#
# Author: M. Josien
# Date:  2021-12-16
#
# Copyright : see License.txt

import sac_de_billes
import merope


NList = [2**i for i in range(4,11)]
voxelRuleList = [merope.VoxelRule.Center, merope.VoxelRule.Average]
voxelRuleNames = ["Center", "Average"]

HGRule = merope.HomogenizationRule
homogRuleList = [HGRule.Largest, HGRule.Smallest, HGRule.Reuss, HGRule.Voigt]
homogRuleNames = ["Largest", "Smallest", "Reuss", "Voigt"]


#--------------------------------------------------------------
# Physics
lambda_gold = 317
lambda_lead = 35.3
lambda_water = 0.606
all_lambdas = [lambda_water, lambda_lead, lambda_gold]

# Dimensions of the box
L0 = 10.
L = [L0, L0, L0]

