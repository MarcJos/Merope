#!/usr/bin/python
# coding: utf-8
# 2024-03
# For building pointers to C function in python
#
# 
# Copyright : see License.txt

from numba import cfunc, types
import merope

class function_py2c:
    """
    Class for automatically build a C function from a python function.
    """
    def __init__(self, lambda_function, dim_from):
        """
        @param lambda_function : a function R^d -> R (if d>1, R^d representer as a list of double)
        @param dim_from : dimension d
        """
        if(dim_from == 1):
            signature = types.float64(types.float64) 
        else:
            signature = types.float64(types.CPointer(types.float64)) # double(double*)
        self.__c_function = cfunc(signature, nopython=True)(lambda_function)
        self.__interf_funcPointer = merope.Interf_FuncPointer(self.__c_function.address, [dim_from, 1])

    def get_funcPointer(self):
        """
        @return a Merope object "Interf_FuncPointer"
        """
        return self.__interf_funcPointer


class texture_py2c:
    """
    Class for automatically build a C function from a python function.
    """
    def __init__(self, lambda_function, dim_from):
        """
        @param lambda_function : a function R^d -> R (if d>1, R^d representer as a list of double)
        @param dim_from : dimension d
        """
        if(dim_from == 1):
            signature = types.float64(types.float64, types.int64) 
        else:
            signature = types.float64(types.CPointer(types.float64), types.int64) # double(double*)
        self.__c_function = cfunc(signature, nopython=True)(lambda_function)
        self.__interf_TexturePointer = merope.Interf_TexturePointer(self.__c_function.address, [dim_from, 1])

    def get_funcPointer(self):
        """
        @return a Merope object "Interf_TexturePointer"
        """
        return self.__interf_TexturePointer


