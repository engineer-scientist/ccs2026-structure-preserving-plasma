cd ~/grokbot-sp/p3; D=$1; shift
[ -f $D/reduced.csv ] || python3 reduce3.py $D > /tmp/gb_red.csv && cat ${D}/reduced.csv 2>/dev/null || cat /tmp/gb_red.csv
for f in $D/detail_*.py; do [ -f "$f" ] && { echo "== $f"; python3 perrank.py $f DMPlexLocatePoints TSStep | grep -v '^   '; }; done
