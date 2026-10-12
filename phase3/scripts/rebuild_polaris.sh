module swap PrgEnv-nvidia PrgEnv-gnu >/dev/null 2>&1; module load cuda/13.0 craype-accel-nvidia80 >/dev/null 2>&1; module unload darshan xalt >/dev/null 2>&1
cd ~/grokbot-sp/ex2gb && rm -f ex2 ex2.o && make PETSC_DIR=$HOME/grokbot-sp/petsc PETSC_ARCH=arch-polaris-cuda-opt ex2 > make.log 2>&1; echo rc=$?
ldd ex2 | grep -iE "mpi|not found"
