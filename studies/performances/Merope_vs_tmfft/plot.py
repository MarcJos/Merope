#!/usr/bin/python
# -*- coding: utf-8 -*- 
# 
# plot.py
# Marc Josien
# 01/10/2021
#
# Compare performances of tmfft and mérope
#

import matplotlib.pyplot as plt
import csv

## extract data
nVox0 = []
nbSpheres = []
timeMerope= []

with open('Time.txt', newline='') as csvfile:
    result_read = csv.reader(csvfile, delimiter=' ', quotechar='|')
    next(result_read)
    for line in result_read:
        nVox0.append(int(line[0]))
        nbSpheres.append(int(line[2]))
        timeMerope.append(float(line[4]))


nVox0_ref =  []
nbSpheres_ref = []
timeTMFFT = []

with open('Time_2021_10_05.txt', newline='') as csvfile:
    result_read = csv.reader(csvfile, delimiter=' ', quotechar='|')
    next(result_read)
    for line in result_read:
        nVox0_ref.append(int(line[0]))
        nbSpheres_ref.append(int(line[2]))
        timeTMFFT.append(float(line[5]))

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
    tTMFFT = []
    for i, nV in enumerate(nVox0):
        if (nV == nVoxRef):
            if(i < len(nbSpheres_ref) and nbSpheres[i] == nbSpheres_ref[i] and nV == nVox0_ref[i]):
                tTMFFT.append(timeTMFFT[i])
                nbSph.append(nbSpheres[i])
                tMerope.append(timeMerope[i])

    plt.plot(nbSph, tMerope, linestyle='-', lw=3, color=colors[j], label = "Merope, nVox = "+str(nVoxRef)+"^3", marker='d', markersize=pointsize)
    plt.plot(nbSph, tTMFFT, linestyle=':', lw=3, color=colors[j], label = "TMFFT, nVox = "+str(nVoxRef)+"^3", marker='d', markersize=pointsize)


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
    tMerope = []
    tTMFFT = []
    speedUpRatio = []
    for i, nV in enumerate(nVox0):
        if (nV == nVoxRef):
            if(i < len(nbSpheres_ref) and nbSpheres[i] == nbSpheres_ref[i] and nV == nVox0_ref[i]):
                tTMFFT.append(timeTMFFT[i])
                nbSph.append(nbSpheres[i])
                tMerope.append(timeMerope[i])
                speedUpRatio.append(timeTMFFT[i]/timeMerope[i])
    plt.plot(nbSph, speedUpRatio, linestyle='-', lw=3, color=colors[j], label = "nVox = "+str(nVoxRef)+"^3", marker='d', markersize=pointsize)


plt.xlabel('Nb Of spheres', fontsize=axislabelsize)
plt.ylabel('Speed-up ratio', fontsize=axislabelsize)

plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=legendlabelsize)
plt.yscale('log')
plt.xscale('log')
plt.savefig('Benchmark_time.pdf')

plt.show()

