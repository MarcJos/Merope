#!/usr/bin/python
# coding: utf-8
# 2022-02
# For creating flux mimicking Dirichlet and Neumann BC
#
# 
# Copyright : see License.txt


import numpy as np
from numpy.core._multiarray_umath import dtype
from enum import Enum

from vtk import vtkStructuredPoints
from vtk import vtkStructuredPointsReader
from vtk import vtkStructuredPointsWriter
from vtk.util import numpy_support as VN

import vtk_manipulations as vtk_manip
import text_manipulation
import merope


class BoundaryConditions(Enum):
    Periodic = 0
    Dirichlet = 1
    Neumann = 2


def prepareGeometry(direction_i, boundaryConditions, vtkFileName, coeffVal_inputFileName="Coeffs.txt",
                    materialId_fileName="material_Sym.vtk", zones_fileName="zone_Sym.vtk"):
    """
    Prepares the geometrical setting (=vtk) for launching amitex-fftp
    :param direction_i: direction in which the corrector is computed
    :param boundaryConditions: boundary conditions to be applied
    :param vtkFileName: initial geometry
    :param coeffVal_inputFileName: values of the (isotropic) coefficient
    :param materialId_fileName: file containing indices denoting each symmetrized zone.
        indice == 0 => flux not reverted
        indice == 1 => flux reverted
    :param zones_fileName: file containing the symmetrized zones
    :param coeffVal_oututFileName: negative flux
    :return: void
    """
    vtkSP_input = vtk_manip.file_to_vtkSP(vtkFileName)
    verifyDimensions(direction_i, boundaryConditions, vtkSP_input)
    nbVox = vtk_manip.dimensions_to_nVox(vtkSP_input.GetDimensions())
    L = vtk_manip.L_from(vtkSP_input.GetSpacing(), nbVox)
    ###
    nbSym, signFlux = signOfFlux(direction_i, boundaryConditions)
    ###
    text_manipulation.multiplyFile(coeffVal_inputFileName, -1)
    text_manipulation.multiplyFile(coeffVal_inputFileName, 0)
    ###
    signFlux = turn_sign_into_material(signFlux) # change signs to materials
    print_dilation(nbVox, L, signFlux, materialId_fileName)
    merope.vox.symmetrize(vtkFileName, zones_fileName, nbSym)


def verifyDimensions(direction_i, boundaryConditions, vtkSP_input):
    DIM = len(boundaryConditions)
    if DIM not in [2,3]:
        raise Exception("Incorrect size of boundary conditions")
    if direction_i not in range(1, 1 + DIM):
        raise Exception("Impossible direction")
    if DIM != len(vtkSP_input.GetDimensions()):
        raise Exception("Incompatible dimensions")


def turn_sign_into_material(my_array):
    def transform_sign(x):
        if x == -1:
            return 1
        else:
            return 0
    return np.vectorize(lambda s :  transform_sign(s))(my_array)


def signOfFlux(direction_i, boundaryConditions):
    DIM = len(boundaryConditions)
    transfomations = [transformation_pair(direction_i, j + 1, boundaryConditions[j]) for j in range(0, DIM)]
    nbSym = [log_2(len(t)) for t in transfomations]
    result_sign = np.zeros(tuple(len(t) for t in transfomations), dtype=int)
    ###
    ijk = [0 for i in range(0,DIM)]
    signs = [0 for i in range(0,DIM)]
    for ijk[0], signs[0] in enumerate(transfomations[0]):
        for ijk[1], signs[1] in enumerate(transfomations[1]):
            if DIM == 2:
                result_sign[ijk[0]][ijk[1]] = signs[0] * signs[1]
            else:
                for ijk[2], signs[2] in enumerate(transfomations[2]):
                    result_sign[ijk[0]][ijk[1]][ijk[2]] = signs[0] * signs[1] * signs[2]
    return nbSym, result_sign


def log_2(n):
    if n == 1:
        return 0
    elif n == 2:
        return 1
    elif n == 4:
        return 2
    else:
        raise Exception("Unexpected")


def transformation_pair(direction_i, direction_j, pair_boundaryConditions):
    if pair_boundaryConditions[0] == BoundaryConditions.Periodic:
        return [transformation(direction_i, direction_j, pair_boundaryConditions[0])]
    else:
        if pair_boundaryConditions[0] == pair_boundaryConditions[1]:
            return [1, transformation(direction_i, direction_j, pair_boundaryConditions[0])]
        else:
            sign0 = transformation(direction_i, direction_j, pair_boundaryConditions[0])
            sign1 = transformation(direction_i, direction_j, pair_boundaryConditions[1])
            return [1, sign1, sign1 * sign0, sign0]


def transformation(direction_i, direction_j, boundarCondition):
    if boundarCondition == BoundaryConditions.Periodic:
        return 1
    if boundarCondition == BoundaryConditions.Dirichlet:
        if direction_i != direction_j:
            return -1
        else:
            return 1
    if boundarCondition == BoundaryConditions.Neumann:
        if direction_i != direction_j:
            return 1
        else:
            return -1


def print_dilation(nbVox, L, valueArray, materialId_fileName):
    """
    Dilate a given (small) array of ints to the scale of a voxellation
    :param nbVox: number of voxels of each copy
    :param valueArray: numpy array containing values between 0 and n for the materialId
    :param materialId_fileName: name of the output vtk
    :return: void
    """
    origin = [0 for l in L]
    numpyDilated = np_dilation(nbVox, valueArray)
    DIM = len(L)
    nbVox_dilated = [nbVox[i] * valueArray.shape[i] for i in range(0, DIM)]
    scalars = vtk_manip.get_scalars_from_numpy(numpyDilated)
    spacing = vtk_manip.spacing_from(L, nbVox)
    #
    vtk_manip.write_vtk_from_vtkSP(
        vtk_manip.get_vtkSP(scalars, nbVox_dilated, origin, spacing), materialId_fileName
    )


def np_dilation(nbVox, valueArray):
    DIM = len(nbVox)
    if  DIM != valueArray.ndim or DIM not in [2,3]:
        raise Exception("The dimensions of the nbVox and valueArray should be compatible")

    myShape = tuple([nbVox[DIM - 1 - i] * valueArray.shape[DIM - 1 - i] for i in range(0, DIM)])
    result = np.zeros(myShape, dtype=np.ushort)
    # indices
    ijk = np.zeros(DIM, dtype=int)
    ijk_min = np.zeros(DIM, dtype=int)
    ijk_max = np.zeros(DIM, dtype=int)

    if DIM == 3:
        for ijk[0] in range(0,valueArray.shape[0]):
            for ijk[1] in range(0, valueArray.shape[1]):
                for ijk[2] in range(0, valueArray.shape[2]):
                    for i in range(0,DIM):
                        ijk_min[i] = nbVox[i] * ijk[i]
                        ijk_max[i] = nbVox[i] * (ijk[i] + 1)
                    result[ijk_min[2]:ijk_max[2], ijk_min[1]:ijk_max[1], ijk_min[0]:ijk_max[0]] \
                        = valueArray[ijk[0]][ijk[1]][ijk[2]]
    if DIM == 2:
        for ijk[0] in range(0, valueArray.shape[0]):
            for ijk[1] in range(0, valueArray.shape[1]):
                    for i in range(0, DIM):
                        ijk_min[i] = nbVox[i] * ijk[i]
                        ijk_max[i] = nbVox[i] * (ijk[i] + 1)
                    result[ijk_min[1]:ijk_max[1], ijk_min[0]:ijk_max[0]] \
                        = valueArray[ijk[0]][ijk[1]]
    # reshape the result
    result.shape = result.size
    return result
