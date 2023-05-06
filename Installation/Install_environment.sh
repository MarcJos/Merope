#!/bin/bash
#
# Environment definition for installation
#
###############################################################################################

if [[ "${HOSTNAME}" == *"pleiades"* ]]
then
    export MEROPE_PYBIND_REPO=https://www-git-cad.intra.cea.fr/DEC/collaboratif/mj263790/copy_pybind 
    export MEROPE_VOROPP_REPO=https://www-git-cad.intra.cea.fr/DEC/collaboratif/mj263790/copy_of_voro-plus-plus

    export MKLROOT="/soft/commun/TOOLS/INTEL/ONEAPI/2021.3.0/buster/mkl/2021.3.0"
    export MKL_ROOT_LIB=${MKLROOT}/lib/intel64
    export MKL_INCLUDE_DIRS=${MKLROOT}/include
    source /soft/pleiades/testing/compilers/gcc/gcc-12.1/buster/x86_64/env.sh 
    source /soft/pleiades/codes/TOOLS/cmake-3.16.1-Linux-x86_64/env.sh
else
    export MEROPE_PYBIND_REPO=https://github.com/pybind/pybind11
    export MEROPE_VOROPP_REPO=https://math.lbl.gov/voro++/download/

    echo "------------------------------"
    echo "------------------------------"
    echo "Please verify the location of MKL"
    echo "------------------------------"
    echo "------------------------------"
    export MKL_ROOT_LIB="/usr/lib/x86_64-linux-gnu"
    export MKL_INCLUDE_DIRS="/usr/include/mkl" 
fi


VOROPP_NAME_DIR=voro-plus-plus




