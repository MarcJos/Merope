#!/usr/bin/python
# coding: utf-8
# 2022-03
# Meant to write xml files for feeding amitex
#
# Copyright : see License.txt
#

def write_xml_header(fileName):
    """
    :return: write the header of a standard amitex input
    """
    with open(fileName, "w") as fic:
        fic.write('''<?xml version="1.0" encoding="UTF-8"?>''' + "\n")
        fic.write('''<!-- Automatically generated by amitex_xml_writer.py . Consider it before modifying it. -->''' + "\n" + "\n")


class Parameters_Fourier_iso:
    """
    Parametrization of the law Fourier_iso in a single material (for computing thermal conductivity)
    """
    def __init__(self, coeff_fileName="coeffs.txt", flux_fileNames=None, numM=1):
        if flux_fileNames is None:
            flux_fileNames = tuple(["Coeffs.txt" for i in range(0, 3)])
        self.numM = numM
        self.coefficient_values_fileName = coeff_fileName
        self.flux_values_fileNames = [flux_fileNames[i] for i in range(0, 3)]

    @staticmethod
    def __write_CoeffK(index, coefficient_fileName, fic):
        fic.write('''<CoeffK Index="''' + str(index) + '''" Type="Constant_Zone" File="''' + coefficient_fileName + '''" Format="ASCII"/>''' + "\n")

    def write_into(self, fileName):
        with open(fileName, "a") as fic:
            fic.write('''<!-- MATERIAL ''' + str(self.numM) + '''-->''' + "\n")
            fic.write('''<Material numM="''' + str(self.numM) + '''" LibK="" LawK="Fourier_iso_polarization">''' + "\n")
            self.__write_CoeffK(1, self.coefficient_values_fileName, fic)
            for i in range(0, 3):
                self.__write_CoeffK(i+2, self.flux_values_fileNames[i], fic)
            fic.write('''</Material>''' + "\n")
            fic.write("\n")


class Material:
    """
    Parametrization for building a file mat_lin_diffusion.xml
    """
    def __init__(self, coeff_K=10, list_of_param_single_mat=None):
        self.law = "Fourier_iso_polarization"
        self.coefficient_reference_Material = coeff_K
        if list_of_param_single_mat is None:
            list_of_param_single_mat = []
        self.__list_of_single_material = list_of_param_single_mat
        for i, param in enumerate(self.__list_of_single_material):
            param.numM = i + 1

    def __write_reference_material(self, fileName):
        with open(fileName, "a") as fic:
            fic.write('''<!-- REFERENCE MATERIAL -->''' + "\n")
            fic.write('''<Reference_MaterialD K0="''' + str(self.coefficient_reference_Material) + '''"/>''' + "\n")
            fic.write("\n")

    def __write_all_materials(self, fileName):
        for param in self.__list_of_single_material:
            param.write_into(fileName)

    def write_into(self, fileName):
        write_xml_header(fileName)
        with open(fileName, "a") as fic:
            fic.write('''<Materials>''' + "\n" + "\n")
        self.__write_reference_material(fileName)
        self.__write_all_materials(fileName)
        with open(fileName, "a") as fic:
            fic.write('''</Materials>''' + "\n" + "\n")


class Loading_diffusion:
    """
    Parametrization for building a file loading_diffusion.xml
    """
    def __init__(self, show_grad=False, show_flux=False, direction_values=None):
        if direction_values is None:
            direction_values = [0., 0., 0.]
        self.show_grad = show_grad
        self.show_flux = show_flux
        self.direction_values = direction_values

    def __write_show_output(self, fileName):
        with open(fileName, "a") as fic:
            fic.write('''<!-- OUPUT QUANTITIES -->''' + "\n")
            fic.write('''<Output>''' + "\n")
            fic.write('''<vtk_FluxDGradD FluxD="''')
            fic.write(str(int(self.show_flux)))
            fic.write('''" GradD="''')
            fic.write(str(int(self.show_grad)))
            fic.write('''"/>''' + "\n")
            fic.write('''</Output>''' + "\n" + "\n")

    @staticmethod
    def __write_single_load(fileName, direction_j, value):
        if direction_j == 1:
            xyz = "x0"
        elif direction_j == 2:
            xyz = "y0"
        elif direction_j == 3:
            xyz = "z0"
        else:
            raise Exception("Unexpected direction")
        with open(fileName, "a") as fic:
            fic.write('''<''' + xyz + " ")
            fic.write('''Driving="GradD" Evolution="Linear" Value="''' + str(value) + '''"/>''' + "\n")

    def __write_loads(self, fileName):
        with open(fileName, "a") as fic:
            fic.write('''<Loading Tag = "1">''' + "\n")
            fic.write('''<Time_Discretization Discretization="User" Nincr="1" />''' + "\n")
            fic.write('''<Time_List>1</Time_List>''' + "\n")
            fic.write('''<Output_vtkList>1</Output_vtkList>''' + "\n")
        for i in range(0, 3):
            self.__write_single_load(fileName, i+1, self.direction_values[i])
        with open(fileName, "a") as fic:
            fic.write('''</Loading>''' + "\n" + "\n")

    def write_into(self, fileName):
        write_xml_header(fileName)
        with open(fileName, "a") as fic:
            fic.write('''<Loading_Output>''' + "\n" + "\n")
        self.__write_show_output(fileName)
        self.__write_loads(fileName)
        with open(fileName, "a") as fic:
            fic.write('''</Loading_Output>''' + "\n" + "\n")


class Algo_diffusion:
    """
    For building the file algo_diffusion.xml
    """
    def __init__(self):
        pass

    @staticmethod
    def write_into(fileName, convergence_criterion_value=1.e-4):
        write_xml_header(fileName)
        with open(fileName, "a") as fic:
            fic.write('''<Algorithm_Parameters>''' + "\n" + "\n")
            ###
            fic.write('''<Algorithm Type="Basic_Scheme">                     <!-- "Default" (Basic_Scheme) ou "Basic_Scheme" -->''' + "\n")
            fic.write('''<Convergence_Criterion Value="''' + str(convergence_criterion_value) + '''"/>        <!-- "Default" (1e-4) ou valeur réelle positive <1e-3 -->''' + "\n")
            fic.write('''<Convergence_Acceleration Value="True"/>        <!-- "True" ou "False" -->''' + "\n")
            fic.write('''<Nitermax Value="3000"/> <!-- Default is 1000 -->''' + "\n")
            fic.write('''</Algorithm>''' + "\n" + "\n")
            ###
            fic.write('''<Diffusion>''' + "\n")
            fic.write('''<Filter Type="Default"/>''' + "\n")
            fic.write('''<Stationary Value="true"/>''' + "\n")
            fic.write('''</Diffusion>''' + "\n")
            fic.write('''</Algorithm_Parameters>''' + "\n" + "\n")


if __name__ == "__main__":
    mat = Material(list_of_param_single_mat=[Parameters_Fourier_iso()])
    mat.write_into("Test_material.xml")
    load = Loading_diffusion(show_grad=True, show_flux=True, direction_values=[1., 0., 0.])
    load.write_into("loading.xml")
    Algo_diffusion.write_into("algo_default_diffusion.xml")