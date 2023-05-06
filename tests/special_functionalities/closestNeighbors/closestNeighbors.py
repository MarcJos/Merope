# -*- coding:utf8 -*-
#
# Author: M. Josien
# Date: 29/08/2022
#
# Copyright : see License.txt
#
# get the closest neighbors of a list of spheres

# Create the list of spheres


### import the library (should accessible through PYTHON_PATH)
import sac_de_billes
import merope

# Create a list of sphere

L = [16., 16., 16.]

typeAlgo = sac_de_billes.TypeAlgo.RSA
nameShape = sac_de_billes.NameShape.Tore
seed = 0
desiredRPhi = [[1., 0.2], [0.5, 0.1]]
tabPhases = [0, 1]
mindist = 0.

theSpheres = sac_de_billes.throwSpheres_3D(typeAlgo, nameShape, L, seed, desiredRPhi, tabPhases, mindist)

for i, sph in enumerate(theSpheres):
    print(str(i) + " : " + str(sph.center) + " , " + str(sph.radius))

# Get the closest neighbors

closestNeighbors = merope.getClosestNeighbors_3D(L, theSpheres)

for (k,v) in closestNeighbors.items():
    print("Sphere " + str(k) + " is close in the sense of the Voronoi diagramm to " + str(v))

# Verification
i_0 = 0
sphere_0 = theSpheres[i_0]

def distance(sph1, sph2):
    res = 0
    for i in range(3):
        res += (sph1.center[i] - sph2.center[i])**2
    res = res**0.5
    res -= sph1.radius + sph2.radius
    return res

list_distance = [(i, distance(sph, sphere_0)) for i, sph in enumerate(theSpheres)] 
list_distance.sort(key = lambda x : x[1])

print("###########")
print([list_distance[i] for i in range(len(closestNeighbors[i_0]))])
print(closestNeighbors[i_0])
