# Phase 4: sourced by PBS scripts. runx NP PPN TAG L NX DT STEPS FIX [extra]  (explicit global NX, for strong scaling)
BASE="$(cat common_opts_v3.txt)"
runx() { NP=$1; P=$2; TAG=$3; L=$4; NX=$5; DT=$6; NS=$7; FIX=$8; shift 8
  FIXOPT=""; [ "$FIX" = 1 ] && FIXOPT="-gb_wrap_hash_factor 2"
  TMAX=$(awk "BEGIN{print $DT*$NS*1.0001}"); T0=$(date +%s.%N)
  $LAUNCH -n $NP $PPNFLAG $P $BINDFLAGS $EXE $BASE -dm_plex_box_upper $L -dm_plex_box_faces $NX -ts_time_step $DT -ts_max_steps $NS -ts_max_time $TMAX $FIXOPT "$@" -log_view :$OUT/logview_$TAG.txt > $OUT/ex2_$TAG.out 2>&1
  RC=$?
  echo "$TAG np=$NP ppn=$P nx=$NX L=$L dt=$DT steps=$NS fix=$FIX particles=$((NX*200)) rc=$RC Np_first=$(grep '^E:' $OUT/ex2_$TAG.out | head -1 | awk '{print $13}') Np_last=$(grep '^E:' $OUT/ex2_$TAG.out | tail -1 | awk '{print $13}') wall_s=$(awk "BEGIN{print $(date +%s.%N) - $T0}") extra=$*" >> $OUT/summary.txt; }
Lof() { awk "BEGIN{printf \"%.15g\", 3.141592653589793*$1}"; }
st() { echo "$(date -u +%FT%TZ) $(hostname -s) job=${PBS_JOBID%%.*} $PBS_JOBNAME $*" >> ~/grokbot-status.txt; }
