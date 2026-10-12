H=$(hostname); mkdir -p ~/grokbot-sp/ex2gb && cd ~/grokbot-sp/ex2gb && case $(hostname) in aurora*) cp ../ex2_326/ex2.c . ;; *) cp ../ex2/ex2.c . ;; esac && patch ex2.c /tmp/gb_patch.diff && rm /tmp/gb_patch.diff
printf 'include ${PETSC_DIR}/lib/petsc/conf/variables\ninclude ${PETSC_DIR}/lib/petsc/conf/rules\n' > makefile
case $H in
 aurora*) A=arch-aurora-cpu-opt; P3EXE=ex2_326;;
 *) module swap PrgEnv-nvidia PrgEnv-gnu >/dev/null 2>&1; module load cuda/13.0 craype-accel-nvidia80 >/dev/null 2>&1; A=arch-polaris-cuda-opt; P3EXE=ex2;;
esac
make PETSC_DIR=$HOME/grokbot-sp/petsc PETSC_ARCH=$A ex2 > make.log 2>&1; echo "make rc=$?"; tail -2 make.log | cut -c1-150
cd ~/grokbot-sp/p3 && cp run_common.sh run_common.sh.v1 && cp /tmp/gb_rc.sh run_common.sh && rm /tmp/gb_rc.sh && rm -f $P3EXE && ln -s ../ex2gb/ex2 $P3EXE && ls -laL $P3EXE && grep FIXOPT= run_common.sh
