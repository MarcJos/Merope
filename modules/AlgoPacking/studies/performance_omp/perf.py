# -*- coding:utf8 -*-
#
# Date: 18/09/2022
#
# Copyright : see License.txt
#
# Test the performance of the openmp parallelization for AlgoWP

import sac_de_billes
import time

nameFile = "Result.txt"


def test(num_threads, L_0):
    t_0 = time.time()
    ##
    L = [L_0, L_0, L_0]
    sac_de_billes.omp_set_num_threads(num_threads)
    theSpheres = sac_de_billes.throwSpheres_3D(sac_de_billes.TypeAlgo.WP, sac_de_billes.NameShape.Tore, L, 0, [[1, 0.5]], [1], 0)
    ##
    t_1 = time.time() - t_0
    with open(nameFile, "a") as fic:
        my_string = str(len(theSpheres)) + " " + str(sac_de_billes.omp_get_num_threads()) + " " + str(t_1) + "\n"
        fic.write(my_string)
        for i in range(10):
            print("##############")
        print(my_string)
        for i in range(10):
            print("##############")


L_0_list = [10, 20, 30, 40, 50, 75, 100]
nbThreads = [1, 2, 4, 8, 16, 32]


for L_0 in L_0_list:
    for num_threads in nbThreads:
        test(num_threads, L_0)
