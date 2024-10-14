# -*- coding:utf8 -*-
#
# Compute time
#
# Author: M. Josien
# Date: 27/09/2021
#

import time

with open("time.txt", "w") as fic:

    tic_0 = time.time()

    import buildVoxellation
    fic.write(str(time.time()-tic_0))
    fic.write("\n")
    tic_0 = time.time()

    import Thermal_amitex
    fic.write(str(time.time()-tic_0))
    fic.write("\n")
    tic_0 = time.time()

    import Thermal_tmfft
    fic.write(str(time.time()-tic_0))
    fic.write("\n")
    tic_0 = time.time()
