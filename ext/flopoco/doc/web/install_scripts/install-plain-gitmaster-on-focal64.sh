#!/bin/bash

yes | sudo apt update && sudo DEBIAN_FRONTEND=noninteractive apt install -y subversion git cmake sollya wget g++ libsollya-dev flex bison libboost-all-dev autotools-dev autoconf automake f2c libblas-dev liblapack-dev libtool liblpsolve55-dev lp-solve

BASEDIR=$PWD
git clone https://github.com/fixif/WCPG.git && cd WCPG && sh autogen.sh && ./configure && make && sudo make install && cd $BASEDIR

#Finally FloPoCo itself, 
git clone https://gitlab.inria.fr/fdupont/flopoco.git

cd flopoco && mkdir build && cd build && cmake .. && make &&  cd $BASEDIR

# build the html documentation in doc/web. 
cd flopoco/build
./flopoco BuildHTMLDoc

# Pure luxury: bash autocompletion. This should be called only once
./flopoco BuildAutocomplete
mkdir ~/.bash_completion.d/
mv flopoco_autocomplete ~/.bash_completion.d/flopoco
echo ". ~/.bash_completion.d/flopoco" >> ~/.bashrc

# Now show the list of operators
./flopoco  
