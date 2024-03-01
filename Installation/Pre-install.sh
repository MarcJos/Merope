#!/bin/bash
#
# Script for downloading the prerequisite
#
# Should be executed in the root of the directory
#
# marc.josien@cea.fr
#
###############################################################################################

source Installation/Install_environment.sh
cd modules

## get  pybind11
cd merope_core/Interface_python
rm -rf pybind11
git clone $MEROPE_PYBIND_REPO pybind11
cd ../../
cp -r merope_core/Interface_python/pybind11 AlgoPacking/Interface_python/
git clone $EIGEN_REPO local_eigen
cp -r local_eigen/Eigen .
rm -rf local_eigen

if [[ "${HOSTNAME}" == *"pleiades"* ]]
then
    ### get voro++
    rm -rf $VOROPP_NAME_DIR
    git clone $MEROPE_VOROPP_REPO $VOROPP_NAME_DIR
    ## put CMake install inside voro++
else
    for i in {1..5}
    do
        echo "------------------------------------------------------------------------------"
    done
    echo " WARNING : "
    echo " Please dowload the sources of voro++ and put them into the folder : modules/"${VOROPP_NAME_DIR}
    echo " The sources of voro++ can be found here : "$MEROPE_VOROPP_REPO 
    echo "then re-execute Pre-install.sh"
    for i in {1..5}
    do
        echo "------------------------------------------------------------------------------"
    done
fi
rm $VOROPP_NAME_DIR/CMakeLists.txt
rm $VOROPP_NAME_DIR/src/CMakeLists.txt
cp cmake/for_voropp/voropp_1.cmake $VOROPP_NAME_DIR/CMakeLists.txt
cp cmake/for_voropp/voropp_2.cmake $VOROPP_NAME_DIR/src/CMakeLists.txt

