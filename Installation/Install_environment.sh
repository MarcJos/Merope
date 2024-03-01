#!/bin/bash
#
# Environment definition for installation
#
###############################################################################################

if [[ "${HOSTNAME}" == *"pleiades"* ]]
then
    export MEROPE_PYBIND_REPO=https://www-git-cad.intra.cea.fr/DEC/collaboratif/mj263790/copy_pybind 
    export MEROPE_VOROPP_REPO=https://www-git-cad.intra.cea.fr/DEC/collaboratif/mj263790/copy_of_voro-plus-plus
    export EIGEN_REPO=https://www-git-cad.intra.cea.fr/DEC/collaboratif/mj263790/eigen_copie

    
    CODENAME=$(lsb_release -c |awk '{print $2}')
    #MKL_VERSION="2021.3.0"
    #export MKLROOT="/soft/commun/TOOLS/INTEL/ONEAPI/${MKL_VERSION}/${CODENAME}/mkl/latest"
    MKL_VERSION="2023.1"
    export MKLROOT="/soft/commun/TOOLS/INTEL/ONEAPI/${MKL_VERSION}/mkl/latest"
    

    source /soft/pleiades/testing/compilers/gcc/gcc-12.1/${CODENAME}/x86_64/env.sh 
    source /soft/pleiades/codes/TOOLS/cmake-3.16.1-Linux-x86_64/env.sh
else
    export MEROPE_PYBIND_REPO=https://github.com/pybind/pybind11
    export MEROPE_VOROPP_REPO=https://math.lbl.gov/voro++/download/
    export EIGEN_REPO=https://gitlab.com/libeigen/eigen

    echo "------------------------------"
    echo "------------------------------"
    echo "Please verify the location of MKL"
    echo "------------------------------"
    echo "------------------------------"
    export MKL_ROOT_LIB="/usr/lib/x86_64-linux-gnu"
    export MKL_INCLUDE_DIRS="/usr/include/mkl" 
fi


VOROPP_NAME_DIR=voro-plus-plus




