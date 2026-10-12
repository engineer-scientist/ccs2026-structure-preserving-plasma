#!/bin/bash
R=/workspace/research/results/pusher
S="-snes_rtol 1e-14 -snes_atol 1e-15 -snes_stol 1e-16"
run(){ f=$1; m=$2; dt=$3; n=$4; u0=$5; extra=$6; ./pusher -field $f -method $m -dt $dt -steps $n -nout $((n/2000>0?n/2000:1)) -u0 $u0 -o $R/${f}_${m}.dat $S $extra > $R/${f}_${m}.log 2>&1; }
UES=1.0,0.5,0,0.5,0,0.3
UUN=1.0,0,0,0,1.0,0.2
USO=1.1,0,0,0,0.004,0.003
for m in boris rk4 dg; do run es $m 0.1 1000000 $UES & done
run es bsi 0.1 1000000 $UES "-ts_basicsymplectic_type 2" &
wait
for m in boris rk4 dg; do run uniform $m 0.2 1000000 $UUN & done
run uniform bsi 0.2 2000 $UUN "-ts_basicsymplectic_type 2" &
for m in boris rk4 dg; do run solovev $m 0.2 1000000 $USO & done
run solovev bsi 0.2 2000 $USO "-ts_basicsymplectic_type 2" &
wait
for f in es uniform solovev; do
  u=$UES; [ $f = uniform ] && u=$UUN; [ $f = solovev ] && u=$USO
  dt=0.2; [ $f = es ] && dt=0.1
  for m in boris rk4 dg bsi; do ./pusher -field $f -method $m -dt $dt -jac 100 -jac_eps 1e-5 -u0 $u $S -ts_basicsymplectic_type 2 2>/dev/null | grep JACDET; done
done > $R/jacdet.txt
echo ALLDONE
