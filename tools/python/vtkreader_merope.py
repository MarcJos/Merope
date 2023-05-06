# -*- coding:utf8 -*-
#
# Exporting vtk file as numpy array, and compare the same field of 2 different vtk files
#
# Date: 07/2021
# Copyright : see License.txt
# 
# inspired by https://stackoverflow.com/questions/11727822/reading-a-vtk-file-with-python

import numpy
import text_manipulation
from vtk import vtkStructuredPointsReader
from vtk.util import numpy_support as VN


def read_file(fileName):
    """
    @param fileName: .vtk file
    @return: a vtkStructuredPoint object
    """
    reader = vtkStructuredPointsReader()
    reader.SetFileName(fileName)
    reader.ReadAllVectorsOn()
    reader.ReadAllScalarsOn()
    reader.Update()
    return reader.GetOutput()


def read_field(fileName, fieldName):
    '''
    @return : numpy field of double of the selected proprerty
    @param fileName : .vtk file
    @param fieldName : name of the selected property
    '''
    all_data = read_file(fileName)
    data = VN.vtk_to_numpy(all_data.GetCellData().GetArray(fieldName))
    dimensions = all_data.GetDimensions()
    ### vtk arrays are stored in the order z,y,x
    dimensions = tuple(reversed([dim-1 for dim in dimensions]))
    data.shape = dimensions
    if len(dimensions) == 2:
        data = numpy.swapaxes(data, 0, 1)
    else:
        data = numpy.swapaxes(data, 0, 2)
    return data


def read_L(fileName):
    reader = vtkStructuredPointsReader()
    reader.SetFileName(fileName)
    reader.Update()
    vtkSP = reader.GetOutput()
    dimensions = vtkSP.GetDimensions()
    spacing = vtkSP.GetSpacing()
    return [spacing[i] * (dimensions[i]-1) for i in range(0, len(spacing))]


def read_nVox(fileName):
    reader = vtkStructuredPointsReader()
    reader.SetFileName(fileName)
    reader.Update()
    vtkSP = reader.GetOutput()
    dimensions = vtkSP.GetDimensions()
    ### vtk arrays are stored in the order z,y,x
    return [dim - 1 for dim in dimensions]


def read_coefficientField(vtk_fileName, coeff_fileName):
    '''
    @param vtk_fileName: vtk file describing each zone
    @param coeff_fileName: .txt file describing for each value of the zone, the value of the coefficient field
    @return: the coefficient field, from the Merope format
    '''
    materials = read_field(vtk_fileName, "MaterialId")
    theShape = materials.shape
    ###
    materials = materials.reshape(-1)
    coefficients = text_manipulation.get_coeff(coeff_fileName)
    result = numpy.array([coefficients[mat] for mat in materials])
    result = result.reshape(theShape)
    ###
    return result


def read_property(fileName, fieldName, renormalize=False, Slice=()):
    '''
    @return : numpy field of ints of the selected proprerty
    @param fileName : .vtk file
    @param fieldName : name of the selected property
    @param Slice : the part of the domain to read
    @param renomalize : boolean, if yes, the minimal phase is always 0 (due to the fact that some vtkwriters use 1 as the minimal phase whereas others use 0)
    '''
    all_data = read_file(fileName)
    data = VN.vtk_to_numpy(all_data.GetCellData().GetArray(fieldName))
    data = numpy.int_(data)

    dimension = all_data.GetDataDimension() # dimensions of the array
    nVox =[n-1 for n in all_data.GetDimensions()]

    if Slice!=():
        Nmin, Nmax = Slice
        slice_data = []
        for y in range(Nmin[1], Nmax[1]):
            for x in range(Nmin[0], Nmax[0]):
                if (dimension==3):
                    for z in range(Nmin[2], Nmax[2]):
                        slice_data.append(data[x+nVox[0]*y+nVox[0]*nVox[1]*z])
                else:
                    slice_data.append(data[x+nVox[0]*y])
        data = slice_data

    if renormalize:
        minimum = min(data)
        data = data - minimum
    return data


def read_double_property(fileName, fieldName):
    '''
    @return : numpy field of double of the selected proprerty
    @param fileName : .vtk file
    @param fieldName : name of the selected property
    '''
    all_data = read_file(fileName)
    data = VN.vtk_to_numpy(all_data.GetCellData().GetArray(fieldName))
    return numpy.double(data)


def compare(fileName1, fileName2, fieldName, threshold):
    '''
    @return : whether the difference between two field OF INT is above a threshold
    @param fileName1, fileName2 : name of the .vtk files to be compared
    @param threshold : value
    '''
    data1 = read_property(fileName1,fieldName, True)
    data2 = read_property(fileName2,fieldName, True)
    absolute_error = max(abs(data1-data2))
    print(absolute_error)
    return threshold > absolute_error 


def compare_double(fileName1, fileName2, fieldName, threshold):
    '''
    @return : whether the difference between two field OF DOUBLE is above a threshold
    @param fileName1, fileName2 : name of the .vtk files to be compared
    @param threshold : value
    '''
    data1 = read_double_property(fileName1,fieldName)
    data2 = read_double_property(fileName2,fieldName)
    absolute_error = max(abs(data1-data2))
    print(absolute_error)
    return threshold > absolute_error 


def compare_error(fileName1, fileName2, fieldName, threshold):
    '''
    throws an error if the difference between two field OF INT is above a threshold
    @see compare
    '''
    if not(compare(fileName1, fileName2, fieldName, threshold)):
        raise NameError("Over the threshold")


def compare_error_double(fileName1, fileName2, fieldName, threshold):
    '''
    throws an error if the difference between two field OF DOUBLE is above a threshold
    @see compare_double
    '''
    if not(compare_double(fileName1, fileName2, fieldName, threshold)):
        raise NameError("Over the threshold")


def differencePercentage(field1, field2, threshold):
    '''
    @return : the percentage of indices i such that (field1[i] - field2[i])> threshold
    @param field1, field2 : lists of numbers
    '''
    numberOfVoxels = len(field1)
    numberOfDifferentVoxels = 0
    if len(field1) != len(field2):
        raise NameError("Incoherent length")
    for i in range(0,len(field1)):
        if abs(field1[i] - field2[i]) > threshold:
            numberOfDifferentVoxels += 1
    return numberOfDifferentVoxels / numberOfVoxels


def compareVoxellationPercentage(fileName1, fileName2, fieldName, Slice1=(), Slice2=()):
    '''
    @return : the percentage of different voxels wrt fieldName of 2 .vtk files 
    @param fileName1, fileName2 : name of the .vtk files to be compared
    @param Slice1, Slice2 : tuple, if empty, the entire domain is compared, else, only a selected slice ([NminX, NminY, NminZ],[NmaxX, NmaxY, NmaxZ]) is taken
    @param fieldName : name of the field to be compared
    '''
    
    field1 = read_property(fileName1, fieldName, True, Slice1)
    field2 = read_property(fileName2, fieldName, True, Slice2)
    return differencePercentage(field1, field2, 0)


def build_mask(lMin, lMax, L, nVox):
    nMin = [int(round(lMin[i] / L[i] * nVox[i])) for i in range(0, len(L))]
    nMax = [int(round(lMax[i] / L[i] * nVox[i])) for i in range(0, len(L))]
    return tuple([slice(nMin[i], nMax[i]) for i in range(0, len(L))])
