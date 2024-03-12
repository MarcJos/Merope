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
## case want OpenMP for FFT
## cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$mon_adresse -DMEROPE_USE_OPENMP_FOR_FFT=True
cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$mon_adresse
### CASE DEBUG : 
###cmake .. -DPython_ADDITIONAL_VERSIONS=3.8 -DCMAKE_BUILD_TYPE=Debug

make -j
make install
