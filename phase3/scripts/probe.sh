cd ~/grokbot-sp; hostname; ls
echo == jobs; qstat -u $USER
for q in debug debug-scaling prod; do echo == $q; qstat -Qf $q 2>/dev/null | grep -iE "max_run|resources_max|resources_min|enabled|started|max_queued"; done
echo == quota; myquota 2>/dev/null; df -h ~ | tail -1
