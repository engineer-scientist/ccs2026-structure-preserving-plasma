#!/bin/bash
# Phase 4 long pusher runs on box: 10x steps (1e7) vs phase-2 1e6
cd /workspace/research/sp/pusher; R=/workspace/research/phase4/box_pusher
S="-snes_rtol 1e-14 -snes_atol 1e-15 -snes_stol 1e-16 -ts_max_snes_failures -1"
run(){ f=$1; m=$2; dt=$3; n=$4; u0=$5; timeout 6h ./pusher -field $f -method $m -dt $dt -steps $n -nout 5000 -u0 $u0 -o $R/${f}_${m}_long.dat $S > $R/${f}_${m}_long.log 2>&1; echo "$f $m rc=$? $(date)" >> $R/done.txt; }
for m in boris rk4 dg; do run uniform $m 0.2 100000000 1.0,0,0,0,1.0,0.2 & run solovev $m 0.2 100000000 1.1,0,0,0,0.004,0.003 & done
wait; echo ALLDONE >> $R/done.txt
