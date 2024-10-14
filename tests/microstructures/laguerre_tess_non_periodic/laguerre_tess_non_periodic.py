
# -*- coding:utf8 -*-
#
# Author: J-M VANSON
# Date: 23/09/2024
#
# Copyright : see License.txt
#
# test periodicity on Laguerre tesselations

import sac_de_billes
import merope


def normeSquare(x, L):
    res = 0
    for i in range(len(x)):
        res += (x[i]-0.5 * L[i])**2
    return res

### Define the spheres

nameShape = sac_de_billes.NameShape.Tore
L = [10, 10, 10]
seed = 0
NbSpheres = 2000
tabRadii = [1 for i in range(NbSpheres)]
tabPhases = [i % 10 for i in range(NbSpheres)]
mindist = 0


theSpheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.BOOL, nameShape, L, seed, tabRadii, tabPhases, mindist)

for i in range(len(theSpheres)):
    theSpheres[i].phase = i % 10

# without periodicity
voro1 = merope.VoroInterface_3D(L,theSpheres, [False, False, False])
voro1.printCustom("%i\n%q\n%P\n%t", "laguerre_tess_custom_no_pbc.dat") # https://math.lbl.gov/voro++/doc/custom.html

# with periodicity
voro2 = merope.VoroInterface_3D(L,theSpheres)
voro2.drawCellsPov("laguerre_tess_gnuplot_pbc.txt")
voro2.drawGnuPlot("laguerre_tess_pov_pbc.pov")
voro2.printCustom("%i\n%q\n%P\n%t", "laguerre_tess_custom_pbc.dat") # https://math.lbl.gov/voro++/doc/custom.html