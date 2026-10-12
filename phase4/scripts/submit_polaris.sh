set -x
mkdir -p ~/grokbot-sp/p4 && cd ~/grokbot-sp/p4 && tar xzf /tmp/gb_p4.tgz && rm /tmp/gb_p4.tgz
ln -sf ../ex2gb/ex2 ex2; ls -la ex2* ; ls -la ../ex2gb/ex2
cp ../p3/reduce3.py . 
qstat -Qf debug debug-scaling  2>&1 | grep -E "Queue:|max_run|max_queued|resources_max.(walltime|nodect)|resources_min|enabled|started"
J1=$(qsub ss_polaris.pbs); echo "SS $J1"
J2=$(qsub lr_polaris.pbs); echo "LR $J2"
echo "$(date -u +%FT%TZ) polaris phase4 SUBMITTED ss=$J1 lr=$J2" >> ~/grokbot-status.txt
sleep 5; qstat -u $USER; qstat -T -u $USER 2>&1 | tail -6; for j in $J1 $J2; do qstat -f $j 2>/dev/null | grep -E "job_state|comment|estimated"; done
