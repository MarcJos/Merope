source /soft/pleiades/testing/compilers/gcc/gcc-12.1/buster/x86_64/env.sh 

mkdir BUILD-DIR INSTALL-DIR
cd INSTALL-DIR
mon_adresse=$PWD
cd ../BUILD-DIR/
cmake ../ -DCMAKE_BUILD_TYPE=Release -DCMAKE_INSTALL_PREFIX=$mon_adresse
make -j
make install
