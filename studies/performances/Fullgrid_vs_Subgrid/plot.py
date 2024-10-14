#!/usr/bin/python
# -*- coding: utf-8 -*- 
# Author: M. Josien
#
# Copyright : see License.txt

import pandas as pd

# Define the function to convert strings to lists of numbers
def convert_to_numbers(x):
    return list(map(float, x.split()))


# Charge the CSV file
# Specify the separator as a space and the converter for the columns concerned
data_Z = pd.read_csv('TimeVox_Z.res', sep='\s+', converters={'col1': convert_to_numbers, 'col2': convert_to_numbers})
data_X = pd.read_csv('TimeVox_X.res', sep='\s+', converters={'col1': convert_to_numbers, 'col2': convert_to_numbers})


import matplotlib.pyplot as plt

pointsize = 10
axislabelsize = 24
axisticksize = 30
legendlabelsize = 20
colors = ['crimson', 'steelblue', 'forestgreen', 'darkorange', 'mediumpurple', 'sienna', 'orchid',
             'grey', 'yellowgreen', 'darkturquoise']

L_all =  [32, 64, 128, 256, 512, 1024]

def print_data(data, name, linestyle):
    for j in range(len(L_all)):
        L = L_all[j]
        print(L)
        data_L = data[data["nbVox0"]==L]
        nbSph = data_L["nSpheres"]
        accel_fact = [list(data_L["Time_3D_full"])[i] / list(data_L["Time_3D_SubGrid"])[i] for i in range(len(data_L["Time_3D_full"]))]
        plt.plot(nbSph, accel_fact, linestyle=linestyle, lw=3, color=colors[j], label = "acceleration factor for slice on " + name + ", nVox = "+str(L)+"^3", marker='d', markersize=pointsize)


print_data(data_X, "X", "-")
print_data(data_Z, "Z", "--")

plt.xlabel('Nb Of spheres', fontsize=axislabelsize)
plt.ylabel('Acceleration factor', fontsize=axislabelsize)



plt.legend(bbox_to_anchor=(1.04,1), loc="upper left", fontsize=legendlabelsize)
plt.yscale('log')
plt.xscale('log')
plt.savefig('Benchmark_accel.pdf')