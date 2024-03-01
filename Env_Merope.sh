#!/bin/bash
# 
# Marc Josien
#
# Environment for Merope librairies

# 1) Source Mérope Paths
# Get the path of this script (assumed to be in the root)
pushd . > /dev/null
SCRIPT_PATH="${BASH_SOURCE[0]}"
if ([ -h "${SCRIPT_PATH}" ]); then
  while([ -h "${SCRIPT_PATH}" ]); do cd `dirname "$SCRIPT_PATH"`;
  SCRIPT_PATH=`readlink "${SCRIPT_PATH}"`; done
fi
cd `dirname ${SCRIPT_PATH}` > /dev/null
SCRIPT_PATH=`pwd`;
popd  > /dev/null

# get the environment compiler of the installation
source ${SCRIPT_PATH}/Installation/Install_environment.sh
# Merope &  RSA_algo
export PYTHONPATH=${SCRIPT_PATH}/INSTALL-DIR/lib:$PYTHONPATH
# Scripts
export PATH=$PATH:${SCRIPT_PATH}/scripts
### tools
export PYTHONPATH=${SCRIPT_PATH}/tools/python/:$PYTHONPATH

# 2) Source AMITEX-FFTP and/or TMFFT
if [[ "${HOSTNAME}" == *"pleiades"* ]]
then
    release=$(lsb_release -c |awk '{print $2}')
    processor=$(uname -m)
    version=$release/$processor
    ### use TMFFT
    export LD_LIBRARY_PATH=/soft/pleiades/testing/BUILDS/PLEIADES-trunk/PREREQUIS/$version/BOOST/lib:$LD_LIBRARY_PATH
##    TMFFT_INSTALL=/soft/pleiades/testing/BUILDS/TMFFT-master/$version
    TMFFT_INSTALL=/soft/pleiades/testing/BUILDS/TMFFT-trunk/$version

    export PYTHONPATH=$TMFFT_INSTALL/lib/python3.7m/site-packages:$PYTHONPATH
    export LD_LIBRARY_PATH=$TMFFT_INSTALL/lib:$LD_LIBRARY_PATH
    export PATH=$PATH:$TMFFT_INSTALL/bin

    ### use amitex_fftp
    AMITEX_FOLDER=/soft/pleiades/testing/BUILDS/AMITEX_FFTP/${version}/amitex_fftp-v8.17.8
else
    # FILL IN
    AMITEX_FOLDER=/usr/lib/amitex_fft-v8.17.8
fi

### use amitex_fftp
export PATH=${AMITEX_FOLDER}/libAmitex/bin:$PATH
export LD_LIBRARY_PATH=${AMITEX_FOLDER}/libAmitex/src/materiauxK:$LD_LIBRARY_PATH
source ${AMITEX_FOLDER}/env_amitex.sh

