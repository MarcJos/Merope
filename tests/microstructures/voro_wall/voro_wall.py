
# -*- coding:utf8 -*-
#
# Author: J-M. VANSON
# Date: 04/10/2024
#
# Copyright : see License.txt
#
# make voronoi cells with cylidrical walls

import sac_de_billes
import merope
import numpy as np
from math import pi, sqrt

# geometrical properties
Rcyl = 4.              # pellet radius (m)
Hcyl = 8              # pellet height (m)
radius = 1.
periodicity = [False, False, False]

L = [2.*Rcyl, 2.*Rcyl, Hcyl]
coo = np.loadtxt('coo_input.in')

theSpheres = []
for i in range(len(coo)):
   sph = sac_de_billes.Sphere_3D()
   sph.center = [coo[i,0],coo[i,1],coo[i,2]]
   sph.radius = radius
   sph.phase = i
   theSpheres.append(sph)


print("Create geometry")
voro = merope.VoroInterface_3D(L,theSpheres, periodicity)
voro.addWallCylinder(Rcyl,Rcyl,0.,0.,0.,Hcyl,Rcyl)
voro.printCustom("%i\n%q\n%P\n%t", "Cells_data.dat") # https://math.lbl.gov/voro++/doc/custom.html
print("Nb. Cells = ", len(theSpheres))

#lag = merope.LaguerreTess_3D(  L, centers)

#axis = merope.geometry.Segment_3D([[Rcyl, Rcyl, 0.], [Rcyl, Rcyl, Hcyl]])
#cylinder = merope.geometry.Cylinder(axis, Rcyl)

#sSphereInc = merope.microInclusion.SphereInclusions(cylinder)

#mi.setInclusions(lag)
#voxellation = merope.Voxellation_3D(mi)
#nbVox = [256 for l in L]
#voxellation.proceed(nbVox)
#fileCoeff = "Coeffs.txt"
#voxellation.printFile(fileVTK, fileCoeff)