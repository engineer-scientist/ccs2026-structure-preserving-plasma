cd ~/grokbot-sp/p3; qstat -u $USER 2>/dev/null | tail -n +6
for d in $(ls -dt runs_* 2>/dev/null | head -3); do echo "== $d"; cat $d/summary.txt; done
