# -*- coding:utf8 -*-
#
# Exploits the output of amitex_fftp
# Date: 27/05/2021
#
# Copyright : see License.txt

import os
import csv
import interface_amitex_fftp.amitex_wrapper as lancement

coeff_amitex_fileName = "thermalCoeff_amitex.txt"


def readThermalFluxGrad(nameFile):
### transforms an output file .std produced by amitex_fftp
### turns it into a pair [flux, gradient]
    with open(nameFile,"r") as inputFic:
        result_read = csv.reader(inputFic, delimiter=' ', quotechar='|')
        for i, line in enumerate(result_read):
            if i==8:
                flux = [float(line[2]), float(line[3]), float(line[4])]
                gradient = [float(line[5]), float(line[6]), float(line[7])]
                ecart_type_flux = [float(line[8]), float(line[9]), float(line[10])]
                ecart_type_grad = [float(line[11]), float(line[12]), float(line[11])]
    return([flux,gradient])


def printThermalCoeff_auxi(resultFiles):
### reads the 9 components of the homogenized thermal conductivity
### after amitex 3 results files
    if len(resultFiles) != 3 :
        raise NameError("Not the expected number of files!")
    table = []
    for i, resFile in enumerate(resultFiles):
        table.append(readThermalFluxGrad(resultFiles[i]))
    homogenized_matrix = [[-table[j][0][i] for j in range(0,3)] for i in range(0, 3)]
    ### Shows the result
    print("----------------")
    print("Homogenized matrix")
    print("----------------")
    for line in homogenized_matrix:
        print(line)
    print("----------------")
    ####
    with open(coeff_amitex_fileName, "w") as fic:
        for line in homogenized_matrix:
            for coeff in line:
                fic.write(str(coeff))
                fic.write(" ")
            fic.write("\n")
    #####
    return homogenized_matrix


def printThermalCoeff(name_folder=None):
### \return : the homogenized thermal conductivity after amitex computed it (see Lancement_AMITEX.py)
### \param name_folder : name of the folder in which compute_thermal_coeff was launched
    resultFiles = [lancement.name_res_file_std(i) for i in range(1,4)]
    return printThermalCoeff_auxi(resultFiles)

