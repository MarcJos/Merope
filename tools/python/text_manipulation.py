#!/usr/bin/python
# coding: utf-8
# 2022-03
# Simple manipulations on csv files
#
# Copyright : see License.txt

import numpy as np


def get_coeff(fileName):
    table = []
    with open(fileName, "r") as inputFile:
        for inputLine in inputFile:
            table.append(float(inputLine))
    return np.array(table)


def max_from_file(fileName):
    coeffs = get_coeff(fileName)
    coeffs = coeffs.reshape(-1)
    return max(coeffs)


def getPrefix(factor_value, given_prefix=None):
    if given_prefix is None:
        prefix = str(abs(factor_value)) + "_"
        if factor_value < 0:
            prefix = "minus" + prefix
    else:
        prefix = given_prefix
    return prefix


def multiplyFile(coeffVal_inputFileName, factor_value, prefix=None):
    prefix = getPrefix(factor_value, prefix)
    coeffVal_oututFileName = prefix + coeffVal_inputFileName
    with open(coeffVal_inputFileName, "r") as inputFile:
        with open(coeffVal_oututFileName, "w") as outputFile:
            for inputLine in inputFile:
                outputFile.write(str(factor_value * float(inputLine)) + "\n")
    return coeffVal_oututFileName
