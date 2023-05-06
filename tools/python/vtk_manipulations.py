# -*- coding:utf8 -*-
#
# VTK manipulation 
# Date: 02/2022
# 
# Copyright : see License.txt

import numpy as np
import vtk as vtk
from vtk import vtkStructuredPoints
from vtk import vtkStructuredPointsReader
from vtk import vtkStructuredPointsWriter
from vtk.util import numpy_support as VN


def from_discreteCoord_to_index(i, nVox):
    """
    :return : the index of the array corresponding to the discrete coordinate
    :param nVox : 3D number of voxels of the array
    :param i : 3D discrete coordinates of the voxel
    """
    return i[0] + nVox[0] * (i[1] + nVox[1] * i[2])


def dimensions_to_nVox(dims):
    return [n-1 for n in dims] 


def nVox_to_dimension(nVox):
    return tuple([n+1 for n in nVox])


def get_scalars_from_numpy(pythonData2):
    scalars = VN.numpy_to_vtk(pythonData2, array_type = vtk.VTK_UNSIGNED_SHORT)
    scalars.SetName("MaterialId")
    return scalars


def file_to_vtkSP(fileName):
    reader = vtkStructuredPointsReader()
    reader.SetFileName(fileName)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    all_data = reader.GetOutput()
    return all_data


def get_vtkSP(scalars, nVox, origin, spacing):
    vtkSPts = vtkStructuredPoints()
    vtkSPts.GetCellData().SetScalars(scalars)
    vtkSPts.SetDimensions(*nVox_to_dimension(nVox))
    vtkSPts.SetOrigin(*origin)
    vtkSPts.SetSpacing(*spacing)
    return vtkSPts


def write_vtk_from_vtkSP(vtkSPts, outputFileName = "NewPhases.vtk"):
    # See https://python.hotexamples.com/fr/examples/vtk/-/vtkStructuredPointsWriter/python-vtkstructuredpointswriter-function-examples.html
    writer = vtkStructuredPointsWriter()
    writer.SetFileName(outputFileName)
    writer.SetInputData(vtkSPts)
    writer.SetFileTypeToBinary()
    writer.Update()
    writer.Write()


def spacing_from(L, nbVox):
    return [L[i]/nbVox[i] for i in range(0, len(nbVox))]


def L_from(spacing, nbVox):
    return [spacing[i] * nbVox[i] for i in range(0, len(nbVox))]
