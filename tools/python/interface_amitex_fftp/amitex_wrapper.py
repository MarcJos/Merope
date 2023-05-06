# -*- coding:utf8 -*-
#
# Launches AMITEX computation
# Date: 27/05/2021
# 
# Copyright : see License.txt
#

import os
import time
import interface_amitex_fftp.amitex_xml_writer as ami_xml
import create_sym
import text_manipulation
import vtkreader_merope


### Names of files/directories
###

AMITEX  = "amitex_fftp" ### put here the amitex executable


def name_main_dir():
    return "amitex_simulation"


def name_res_dir(i):
    ### \return : name of the result directory for amitex simulation
    ### \param i : spatial direction, i \in {1,2,3}
    return name_main_dir() + "/" + "result_" + str(i)


def name_res_file(i):
    ### \return : name of the result file for amitex simulation
    ### \param i : spatial direction, i \in {1,2,3}
    return name_res_dir(i) + "/" + "output_file"


def name_res_file_std(i):
    return name_res_file(i) + ".std"


def make_res_dir():
    os.system("rm -rf " + name_main_dir())
    os.system("mkdir " + name_main_dir())
    for i in range (1, 4):
        os.system("mkdir " + name_res_dir(i))


def name_flux_file(i, j):
    return name_res_dir(i) + "/output_file_1_" + name_flux_variable(i, j) + ".vtk"

def name_grad_file(i, j):
    return name_res_dir(i) + "/output_file_1_" + name_grad_variable(i, j) + ".vtk"

def name_flux_variable(i, j):
    return "Flux1_" + str(j)

def name_grad_variable(i, j):
    return "gradQ1_" + str(j)

#####
#####


def set_direction(i, coeff_fileName, reflexion_signature=1):
    null_coeff_fileName = text_manipulation.multiplyFile(coeff_fileName, 0)
    activatedFlux_fileName = text_manipulation.multiplyFile(coeff_fileName, -1 * reflexion_signature)
    flux_fileNames = [null_coeff_fileName for i in range(0, 3)]
    flux_fileNames[i-1] = activatedFlux_fileName
    return flux_fileNames


def compute_single_thermal_coeff(direction_i, boundaryConditions=None, zone_vtk="Zones.vtk",
                                 coeff_fileName="Coeffs.txt", nb_proc=2, display_field=False,
                                 convergence_criterion_value=1.0e-4):
    if boundaryConditions is None:
        boundaryConditions = [[create_sym.BoundaryConditions.Periodic] for i in range(0,3)]
    # names of files
    mat_sym_vtk = "material_Sym.vtk"
    zone_sym_vtk = "zone_Sym.vtk"
    mat_xml = "mat_lin_diffusion.xml"
    load_xml = "loading.xml"
    algo_xml = "algo_default_diffusion.xml"
    # algo_xml
    ami_xml.Algo_diffusion.write_into(algo_xml, convergence_criterion_value)
    # maximal coefficient
    coeff_K = text_manipulation.max_from_file(coeff_fileName)
    # loading_xml
    load_conf = ami_xml.Loading_diffusion(show_flux=display_field, show_grad=display_field)
    load_conf.write_into(load_xml)
    # compute
    all_periodic = all([bc[0] == create_sym.BoundaryConditions.Periodic for bc in boundaryConditions])
    if all_periodic:
        flux_fileNames = set_direction(direction_i, coeff_fileName)
        # mat_xml
        param_Fourier = ami_xml.Parameters_Fourier_iso(coeff_fileName=coeff_fileName,
                                                       flux_fileNames=flux_fileNames)
        mat = ami_xml.Material(coeff_K=coeff_K, list_of_param_single_mat=[param_Fourier])
        mat.write_into(mat_xml)
        #
        exec_amitex(zone_vtk, mat_xml, algo_xml, load_xml, direction_i, nb_proc)
    else:
        create_sym.prepareGeometry(direction_i, boundaryConditions, zone_vtk, coeff_fileName, mat_sym_vtk,
                                   zone_sym_vtk)
        if max(vtkreader_merope.read_property(mat_sym_vtk, "MaterialId")) > 0: ### non trivial material
            flux_fileNames = [set_direction(direction_i, coeff_fileName, 1),
                              set_direction(direction_i, coeff_fileName, -1)]
        else:
            flux_fileNames = [set_direction(direction_i, coeff_fileName, 1)]
        # mat_xml
        params_Fourier = [ami_xml.Parameters_Fourier_iso(coeff_fileName=coeff_fileName, flux_fileNames=ffn)
                          for ffn in flux_fileNames]
        mat = ami_xml.Material(coeff_K=coeff_K, list_of_param_single_mat=params_Fourier)
        mat.write_into(mat_xml)
        #
        exec_amitex(zone_sym_vtk, mat_xml, algo_xml, load_xml, direction_i, nb_proc, mat_vtk=mat_sym_vtk)


def computeThermalCoeff(zone_vtk, nb_proc, boundaryConditions=None, coeff_fileName="Coeffs.txt", display_field=False,
                        convergence_criterion_value=1.0e-4):
    # Computes the diffusion
    # Requires 2 things
    # :param zone_vtk : a .vtk file, that contains the description of the zones
    # :param nb_proc : number of processors
    make_res_dir()
    for i in range(1, 4):
        compute_single_thermal_coeff(i, boundaryConditions=boundaryConditions, zone_vtk=zone_vtk, nb_proc=nb_proc,
                                     coeff_fileName=coeff_fileName, display_field=display_field,
                                     convergence_criterion_value=convergence_criterion_value)
   

def exec_amitex(zone_vtk, mat_xml, algo_xml, load_xml, index_result=0, nb_proc=2, mat_vtk=None):
    name_result_file = name_res_file(index_result)
    tic0 = time.time()
    command_amitex = "mpirun -np " + str(nb_proc) + " " + AMITEX + " -nz " + zone_vtk + " -m " + mat_xml + " -a " + algo_xml + " -c " + load_xml + " -s " + name_result_file
    if mat_vtk is not None:
        command_amitex = command_amitex + " -nm " + mat_vtk
    print(command_amitex)
    os.system(command_amitex)
    os.system("""sed -i "s/\ \ /\ /g" """ + name_result_file + ".std")
    print(time.time() - tic0)

