# sourced by phase-3 PBS scripts. run NP PPN TAG L DT STEPS FIX [extra opts]
# FIX=1 adds -gb_wrap_hash_factor 2 (Grok Bot ex2 patch: ~cell-size hash bins on the periodic-wrap rank only)
BASE="$(cat common_opts_v3.txt)"
run() { NP=$1; P=$2; TAG=$3; L=$4; DT=$5; NS=$6; FIX=$7; shift 7; NX=$((256*NP))
  FIXOPT=""; [ "$FIX" = 1 ] && FIXOPT="-gb_wrap_hash_factor 2"
  TMAX=$(awk "BEGIN{print $DT*$NS*1.0001}")
  T0=$(date +%s.%N)
  $LAUNCH -n $NP $PPNFLAG $P $BINDFLAGS $EXE $BASE -dm_plex_box_upper $L -dm_plex_box_faces $NX -ts_time_step $DT -ts_max_steps $NS -ts_max_time $TMAX $FIXOPT "$@" -log_view :$OUT/logview_$TAG.txt${DETAIL:+,:$OUT/detail_$TAG.py:ascii_info_detail} > $OUT/ex2_$TAG.out 2>&1
  RC=$?
  echo "$TAG np=$NP ppn=$P nx=$NX L=$L dt=$DT steps=$NS fix=$FIX particles=$((NX*200)) rc=$RC Np_first=$(grep '^E:' $OUT/ex2_$TAG.out | head -1 | awk '{print $13}') Np_last=$(grep '^E:' $OUT/ex2_$TAG.out | tail -1 | awk '{print $13}') wall_s=$(awk "BEGIN{print $(date +%s.%N) - $T0}") extra=$*" >> $OUT/summary.txt; }
Lof() { awk "BEGIN{printf \"%.15g\", 3.141592653589793*$1}"; }   # grow-the-box: L = pi*P (4 pi * P/4; width pi per rank)
