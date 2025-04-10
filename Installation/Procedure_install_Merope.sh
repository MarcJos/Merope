#!/bin/bash
#
# Script for installing merope
#
# Should be executed in the root of the directory
#
# marc.josien@cea.fr
#
###############################################################################################

### Environment-dependent prequisite
source Installation/Install_environment.sh

mkdir BUILD-DIR INSTALL-DIR
cd INSTALL-DIR
mon_adresse=$PWD
cd ../BUILD-DIR/
MEROPE_BUILD_TYPE=Release
cmake ../ -DCMAKE_BUILD_TYPE=$MEROPE_BUILD_TYPE -DCMAKE_INSTALL_PREFIX=$mon_adresse -DMEROPE_USE_OPENMP_FOR_FFT=True
# -DMEROPE_COVERAGE=OFF

make -j
make install
