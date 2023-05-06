# -*- coding:utf8 -*-
#
# Example for the use of composite voxels.
# We want to compute the conductivity of lead spheres coated with gold inside water.
#
# Author: M. Josien
# Date:  2021-12-16
#

from parametrization import *
import matplotlib.pyplot as plt
import csv

import sac_de_billes
import merope


## plot
pointsize = 10
axislabelsize = 24
axisticksize = 30
legendlabelsize = 20
colors = ['crimson', 'steelblue', 'forestgreen', 'darkorange', 'mediumpurple', 'sienna', 'orchid',
             'grey', 'yellowgreen', 'darkturquoise']
#


listTypes = []
listCoeffs = []

for i in range(0, len(voxelRuleList)):
    voxRule = voxelRuleList[i]
    voxName = voxelRuleNames[i]
    for j in range(0, len(homogRuleList)):
        hgRule = homogRuleList[j]
        hgName = homogRuleNames[j]
        if(i == 1 or j == 0): ## for center, only 1 name
            listTypes.append((voxName,hgName))
            listCoeffs.append({})
            fileName = "Result" + voxName + "_" + hgName + ".res"
            with open(fileName, newline="") as csvfile:
                result_read = csv.reader(csvfile, delimiter=' ', quotechar='|')
                for line in result_read:
                    listCoeffs[-1][int(line[0])] = float(line[1])
print(listTypes)
print(listCoeffs)

for i, typeVox in enumerate(listTypes):
    n = [nvox for nvox in listCoeffs[i].keys()]
    dx = [L0/nvox for nvox in n]
    value = [listCoeffs[i][nvox] for nvox in n]
    plt.plot(dx,value,color=colors[i],label = typeVox)

plt.xlabel('Voxel size', fontsize = legendlabelsize)
plt.ylabel('Computed effective conductivity', fontsize = legendlabelsize)
plt.xticks(fontsize = legendlabelsize)
plt.yticks(fontsize = legendlabelsize)
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=legendlabelsize)

plt.show()


#####################################################
### IMPORTANT : DEPEND ON THE EXPERIMENT
### get reference value
referenceValue = 0

for i, typeVox in enumerate(listTypes):
    if typeVox[0] == voxelRuleNames[1] and typeVox[1] == homogRuleNames[2]:
        nmax = max(listCoeffs[i].keys())
        referenceValue = listCoeffs[i][nmax]
        del listCoeffs[i][nmax]
#####################################################

for i, typeVox in enumerate(listTypes):
    n = [nvox for nvox in listCoeffs[i].keys()]
    dx = [L0/nvox for nvox in n]
    erreur_relative = [abs((listCoeffs[i][nvox] - referenceValue)/referenceValue) for nvox in n]
    plt.plot(dx,erreur_relative,color=colors[i],label = typeVox)

plt.xlabel('Voxel size', fontsize = legendlabelsize)
plt.ylabel('Relative error on the effective conductivity', fontsize = legendlabelsize)
plt.xticks(fontsize = legendlabelsize)
plt.yticks(fontsize = legendlabelsize)
plt.xscale('log')
plt.yscale('log')
plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=legendlabelsize)

plt.show()



