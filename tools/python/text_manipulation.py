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


def check_lines_equal(file1, file2, start_line, end_line):
    """
    Check if the lines from start_line to end_line in file1 and file2 are equal.

    Parameters:
    file1 (str): Path to the first file.
    file2 (str): Path to the second file.
    start_line (int): Starting line number (1-based index).
    end_line (int): Ending line number (1-based index).

    Returns:
    bool: True if the lines are equal, False otherwise.
    """
    try:
        with open(file1, 'r') as f1, open(file2, 'r') as f2:
            lines1 = f1.readlines()
            lines2 = f2.readlines()

            # Ensure the files have enough lines
            if len(lines1) < end_line or len(lines2) < end_line:
                print("One or both files do not have enough lines.")
                return False

            # Extract the specified lines
            lines1_subset = lines1[start_line - 1:end_line]
            lines2_subset = lines2[start_line - 1:end_line]

            # Check if the lines are equal
            return lines1_subset == lines2_subset
    except FileNotFoundError:
        print("One or both files not found.")
        return False

def check_lines_equal_error(file1, file2, start_line, end_line):
    res = check_lines_equal(file1, file2, start_line, end_line)
    if res:
        print("Comparison success : " + file1 + " / " + file2)
    else:
        raise Exception("Different files : " + file1 + " / " + file2)