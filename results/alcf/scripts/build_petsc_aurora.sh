#!/bin/bash -l
# CPU-only optimized PETSc release build for ex2 diagnosis/pilots on Aurora (HOME dir). Grok Bot.
set -e
cd ~/grokbot-sp
[ -d petsc ] || git clone -q --depth 1 -b release https://gitlab.com/petsc/petsc.git petsc
cd petsc
export PETSC_DIR=$PWD PETSC_ARCH=arch-aurora-cpu-opt
./configure --with-cc=mpicc --with-cxx=mpicxx --with-fc=0 --with-debugging=0 COPTFLAGS="-O3" CXXOPTFLAGS="-O3" --with-blaslapack-dir=$MKLROOT --with-x=0 > ../configure_aurora.log 2>&1
make -j16 > ../make_aurora.log 2>&1
echo BUILD_OK >> ../make_aurora.log
