# -*- coding:utf8 -*-
#
# Computing homogenized coefficients by windowing
#
# Date: 03/2022

import vtkreader_merope as vtk_rM
import interface_amitex_fftp.amitex_wrapper as amitex
import numpy as np


def get_windowed_homogCoeff(percentage=0.5, vtkZone_fileName="Phases.vtk"):
    """
    :param percentage: percentage of the volume used while windowing
    :param vtkZone_fileName: name of the vtk file containing the phases (without symmetries)
    :return: homogenized_matrix
    """
    nVox = vtk_rM.read_nVox(vtkZone_fileName)
    DIM = len(nVox)
    ###
    pure_homogCoeff = np.zeros((DIM, DIM))
    average_nabla_T = np.zeros((DIM, DIM))
    for i in range(0, DIM):
        for j in range(0, DIM):
            pure_homogCoeff[i][j] = get_windowed_coordHomogCoeff(i+1, j+1, percentage, vtkZone_fileName)
            average_nabla_T[i][j] = get_windowed_average_gradiant(i+1, j+1, percentage, vtkZone_fileName)
            if (i==j):
                average_nabla_T[i][j] += 1
    homogCoeff = pure_homogCoeff * np.linalg.inv(average_nabla_T)
    for i in range(0,10):
        print("*****************")
    print(pure_homogCoeff)
    for i in range(0,10):
        print("*****************")
    print(homogCoeff)
    for i in range(0,10):
        print("*****************")
    print(average_nabla_T)
    for i in range(0,10):
        print("*****************")
    return homogCoeff


def get_windowed_coordHomogCoeff(i, j, percentage, vtkZone_fileName):
    """
    :param i, j: coordinates of the homogenized matrix
    :param percentage: percentage of the volume used while windowing
    :param vtkZone_fileName: name of the vtk file containing the phases (without symmetries)
    :return: average of flux \int_{Q_l} e_i . a(e_j + \nabla \phi_j)
    """
    flux1_fileName = amitex.name_flux_file(j, i)
    flux_name = amitex.name_flux_variable(j, i)
    #
    L = vtk_rM.read_L(vtkZone_fileName)
    nVox = vtk_rM.read_nVox(vtkZone_fileName)
    DIM = len(nVox)
    #
    flux_corrector = vtk_rM.read_field(flux1_fileName, flux_name)
    mask = get_mask(DIM, percentage, L, nVox)
    return get_homogCoeff(flux_corrector, mask)

def get_windowed_average_gradiant(i, j, percentage, vtkZone_fileName):
    """
    :param i, j: coordinates
    :param percentage: percentage of the volume used while windowing
    :param vtkZone_fileName: name of the vtk file containing the phases (without symmetries)
    :return: average of gradient \int_{Q_l} e_i . (e_j + \nabla \phi_j)
    """
    grad_fileName = amitex.name_grad_file(j, i)
    grad_name = amitex.name_grad_variable(j, i)
    #
    L = vtk_rM.read_L(vtkZone_fileName)
    nVox = vtk_rM.read_nVox(vtkZone_fileName)
    DIM = len(nVox)
    #
    grad_corrector = vtk_rM.read_field(grad_fileName, grad_name)
    mask = get_mask(DIM, percentage, L, nVox)
    return get_grad(grad_corrector, mask)


def get_mask(DIM, percentage, L, nVox):
    lMin = [0 for i in range(DIM)]
    lMax = [0 for i in range(DIM)]
    if isinstance(percentage, float) or isinstance(percentage, int):
        min_L = (1 - percentage ** (1 / DIM)) / 2
        lMin = [min_L * L[i] for i in range(0, DIM)]
        lMax = [L[i] * (1 - min_L) for i in range(0, DIM)]
    elif isinstance(percentage, list):
        for i in range(DIM):
            lMin_i = 0
            lMax_i = L[i]
            if percentage[i][0] > 0:
                lMin_i = percentage[i][0] * L[i]
            if percentage[i][1] < 1:
                lMax_i = percentage[i][1] * L[i]
            lMin[i] = lMin_i
            lMax[i] = lMax_i
    else:
        raise Error("Unexpected")
    return vtk_rM.build_mask(lMin, lMax, L, nVox)

def get_homogCoeff(flux, mask):
    size_of_mask = len(flux[mask].reshape(-1))
    flux_extracted = -flux[mask]
    result = sum(flux_extracted.reshape(-1)) / size_of_mask
    return result


def get_grad(grad, mask):
    size_of_mask = len(grad[mask].reshape(-1))
    grad_extracted = -grad[mask]
    result = sum(grad_extracted.reshape(-1)) / size_of_mask
    return result

if __name__ == "__main__":
    percentage = 0.5
    print(get_windowed_homogCoeff(percentage, vtkZone_fileName="Composite.vtk"))
