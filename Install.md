# Install

See [cea_install.md](https://www-git-cad.intra.cea.fr/DEC/pleiades/merope/merope_nucleaire/-/blob/master/doc/CEA_install.md) for installation in the CEA.

## Stand-alone installation of sac_de_billes

The module sac_de_billes can be installed alone. See folder [AlgoPacking](modules/AlgoPacking) and documentation therein.

## OS

`Mérope` can be compiled on on Debian Buster, Debian Bullseye, Ubuntu Jammy.

## Prerequisites

Mérope needs :
- a recent Linux distribution (Debian 10+, Ubuntu 22+)
- python 3,
- C++ 17,
- cmake 3.16+,
- the MKL library,
- gcc 10.2+ with g++,
- openmp,
- git.

For compiling, you should have access to repositories containing `pybind11` and `voro++`, and dowload them. Please modify if necessary [Install_environment.sh](Installation/Install_environment.sh).

Moreover, for running, Mérope can be used in combination with the following applications :
- TMFFT (CEA only),
- AMITEX_FFTP,
- MFront.

## Installation

- get the Mérope repository :  
`git clone https://www-git-cad.intra.cea.fr/DEC/pleiades/merope/merope_bibliotheque`

### At the CEA
- get associated projects/prerequisites :  
`bash Installation/Pre-install.sh`  
:warning: all the loadings should be done successfully (**no password error**), otherwise, it should be redone.
- install :  
`bash Installation/Procedure_install_Merope.sh`

### From outside the CEA
- get associated projects/prerequisites :  
    `bash Installation/Pre-install.sh`   
    Please follow the instructions displayed in your terminal.
    - pybind is automatically downloaded from [here](https://github.com/pybind/pybind11)
    - voro++ should be manually dowloaded from [here](https://math.lbl.gov/voro++/download/)
- install :  
`bash Installation/Procedure_install_Merope.sh`

See closed issues for various problems that may arise.


## Use
- source the Mérope environment :  
    - `source Env_Merope.sh`  
    :warning: If you do not have access to the Pleiades machines, you should modify `Env_Merope.sh` to define the path to your current installation of AMITEX-FFTP.
- launch the tests :   
    - `source Env_Merope.sh`  
    - `cd tests/microstructures`  
    - `python3 non-regression_tests.py` [Beware, tests may fail depending on the version of gcc.]
:construction:

## Remark

:warning: OneMKL library could be replaced by FFTW when compiling and linking. [https://www.smcm.iqfr.csic.es/docs/intel/mkl/mkl_manual/appendices/mkl_appG_FFTW3_Intro.htm](https://www.smcm.iqfr.csic.es/docs/intel/mkl/mkl_manual/appendices/mkl_appG_FFTW3_Intro.htm)
Nevertheless, please consider that FFTW is under GPL licence, which is not compatible with the licence of Mérope.
