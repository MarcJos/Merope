#!/usr/bin/python
# -*- coding: utf-8 -*- 
# 
# plot.py
# Marc Josien
# 01/10/2021
#
# Copyright : see License.txt
#

# Compare performances of tmfft and mérope
#

import matplotlib.pyplot as plt
import csv

## extract data
nVox0 = []
nbSpheres = []
timeMerope= []
timeNeper = []



with open('TimeVox.res', newline='') as csvfile:
    result_read = csv.reader(csvfile, delimiter=' ', quotechar='|')
    next(result_read)
    for line in result_read:
        nVox0.append(float(line[1]))
        nbSpheres.append(float(line[0]))
        timeMerope.append(float(line[3]))
        timeNeper.append(float(line[2]))


## plot
pointsize = 10
axislabelsize = 24
axisticksize = 30
legendlabelsize = 20
colors = ['crimson', 'steelblue', 'forestgreen', 'darkorange', 'mediumpurple', 'sienna', 'orchid',
             'grey', 'yellowgreen', 'darkturquoise']



#fixme
reference_nVox = [32, 64, 128, 256, 512]

for j, nVoxRef in enumerate(reference_nVox):
    nbSph = []
    tMerope = []
    tNeper = []
    for i, nV in enumerate(nVox0):
        if (nV == nVoxRef):
            nbSph.append(nbSpheres[i])
            tMerope.append(timeMerope[i])
            tNeper.append(timeNeper[i])
    plt.plot(nbSph, tMerope, linestyle='-', lw=3, color=colors[j], label = "Merope, nVox = "+str(nVoxRef)+"^3", marker='d', markersize=pointsize)
    plt.plot(nbSph, tNeper, linestyle=':', lw=3, color=colors[j], label = "Neper, nVox = "+str(nVoxRef)+"^3", marker='d', markersize=pointsize)


plt.xlabel('Nb Of spheres', fontsize=axislabelsize)
plt.ylabel('Execution time', fontsize=axislabelsize)

plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=legendlabelsize)
plt.yscale('log')
plt.xscale('log')
plt.savefig('Benchmark_time.pdf')

plt.show()


reference_nVox = [32, 64, 128, 256, 512]

for j, nVoxRef in enumerate(reference_nVox):
    nbSph = []
    speedUpRatio = []
    for i, nV in enumerate(nVox0):
        if ( abs(nV -nVoxRef) < 0.1):
            nbSph.append(nbSpheres[i])
            speedUpRatio.append(timeNeper[i]/timeMerope[i])
    plt.plot(nbSph, speedUpRatio, linestyle='-', lw=3, color=colors[j], label = "nVox = "+str(nVoxRef)+"^3", marker='d', markersize=pointsize)
    print(speedUpRatio)

plt.xlabel('Nb Of spheres', fontsize=axislabelsize)
plt.ylabel('Speed-up ratio', fontsize=axislabelsize)

plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=legendlabelsize)
plt.yscale('log')
plt.xscale('log')
plt.savefig('Benchmark_time.pdf')

plt.show()
