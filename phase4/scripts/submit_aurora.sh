set -x
mkdir -p ~/grokbot-sp/p4 && cd ~/grokbot-sp/p4 && tar xzf /tmp/gb_p4.tgz && rm /tmp/gb_p4.tgz
ln -sf ../ex2gb/ex2 ex2_326; ls -la ex2* ; ls -la ../ex2gb/ex2
cp ../p3/reduce3.py . 
qstat -Qf debug debug-scaling prod 2>&1 | grep -E "Queue:|max_run|max_queued|resources_max.(walltime|nodect)|resources_min|enabled|started"
J1=$(qsub ss_aurora.pbs); echo "SS $J1"
J2=$(qsub lr_aurora.pbs); echo "LR $J2"
J3=$(qsub ws_aurora_32.pbs); echo "WS32 $J3"
J4=$(qsub ws_aurora_64.pbs); echo "WS64 $J4"
echo "$(date -u +%FT%TZ) aurora phase4 SUBMITTED ss=$J1 lr=$J2 ws32=$J3 ws64=$J4" >> ~/grokbot-status.txt
sleep 5; qstat -u $USER; qstat -T -u $USER 2>&1 | tail -8; for j in $J1 $J2 $J3 $J4; do qstat -f $j 2>/dev/null | grep -E "job_state|comment|estimated"; done
