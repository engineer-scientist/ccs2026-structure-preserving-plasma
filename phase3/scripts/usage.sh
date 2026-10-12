date -u; qstat -u $USER | tail -n +6
for j in "$@"; do echo "== $j"; qstat -xf $j 2>/dev/null | grep -E "job_state|resources_used.walltime|Exit_status|Resource_List.nodect|stime|obittime"; done
