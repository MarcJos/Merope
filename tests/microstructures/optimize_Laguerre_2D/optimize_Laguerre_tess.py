
# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 19/12/2023
#
# Copyright : see License.txt
#
# Optimize a polycrystal

import sac_de_billes
import merope


def plot_tessellation(L, centers, fileVTK):
    lag = merope.LaguerreTess_2D(  L, centers)
    mi = merope.MultiInclusions_2D()
    mi.setInclusions(lag)
    voxellation = merope.Voxellation_2D(mi)
    nbVox = [1024 for l in L]
    voxellation.proceed(nbVox)
    fileCoeff = "Coeffs.txt"
    voxellation.printFile(fileVTK, fileCoeff)

def normeSquare(x, L):
    res = 0
    for i in range(len(x)):
        res += (x[i]-0.5 * L[i])**2
    return res

### Define the spheres

nameShape = sac_de_billes.NameShape.Tore
L = [10, 10]
seed = 0
NbSpheres = 100
tabRadii = [1 for i in range(NbSpheres)]
tabPhases = [i % 10 for i in range(NbSpheres)]
mindist = 0

theSpheres = sac_de_billes.throwSpheres_2D(sac_de_billes.TypeAlgo.BOOL, nameShape, L, seed, tabRadii, tabPhases, mindist)

for i in range(len(theSpheres)):
    theSpheres[i].phase = i % 10

plot_tessellation(L, theSpheres, "original_tessellation.vtk")

### Optimize the volumes with all equal volumes
total_volume = 1
for i in range(len(L)):
    total_volume *= L[i]
desiredVolumes = [total_volume / NbSpheres for i in range(NbSpheres)]

algo = merope.algo_fit_volumes_2D(L, theSpheres, desiredVolumes)
max_delta_Volume = 1e-6 * total_volume
max_iter = 3000
verbose = True
print("The maximal error on the crystallite volumes before optimization is " + str(algo.maxDeltaVolumes()))
algo.proceed(max_delta_Volume, max_iter, verbose)
print("The maximal error on the crystallite volumes after optimization is " + str(algo.maxDeltaVolumes()))

the_new_spheres = algo.getCenterTessels()
plot_tessellation(L, the_new_spheres, "all_tessels_with_same_volume.vtk")

### Optimize the volumes with different volumes
weighting_function = lambda x : normeSquare(x,L)
desiredVolumes = [weighting_function(sphere.center) for sphere in theSpheres]
total_weight = 0
for weight in desiredVolumes:
    total_weight += weight

desiredVolumes = [weight / total_weight * total_volume for weight in desiredVolumes]

algo = merope.algo_fit_volumes_2D(L, theSpheres, desiredVolumes)
max_delta_Volume = 1e-6 * total_volume
max_iter = 3000
verbose = True
print("The maximal error on the crystallite volumes before optimization is " + str(algo.maxDeltaVolumes()))
algo.proceed(max_delta_Volume, max_iter, verbose)
print("The maximal error on the crystallite volumes after optimization is " + str(algo.maxDeltaVolumes()))

the_new_spheres = algo.getCenterTessels()
plot_tessellation(L, the_new_spheres, "tessels_with_weight.vtk")