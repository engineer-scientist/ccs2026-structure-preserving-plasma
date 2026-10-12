#!/bin/bash
# run on ALCF login: bash ~/grokbot-sp/p4/collect.sh -> ~/grokbot-sp/p4/p4_results_<host>.tgz (no big E: logs >50MB)
cd ~/grokbot-sp/p4 && H=$(hostname -s | cut -d- -f1)
tar czf p4_results_$H.tgz --exclude='*.tgz' runs_* gb_p4_*.o* *.pbs 2>/dev/null; ls -la p4_results_$H.tgz
grep -h . runs_*/summary.txt 2>/dev/null; tail -20 ~/grokbot-status.txt
