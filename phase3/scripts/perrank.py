# Parse PETSc ascii_info_detail python log; print per-rank time/flop for selected events
import sys, re
EV = sys.argv[2:] or ["DMPlexLocatePoints","DMSwarmMigrate","TSStep","SNESSolve"]
pat = re.compile(r'Stages\["Main Stage"\]\["(\w+)"\]\[(\d+)\] = \{"count" : (\d+), "time" : ([0-9.e+-]+).*"flop" : ([0-9.e+-]+)\}')
d = {}
for line in open(sys.argv[1]):
    m = pat.match(line)
    if m and m.group(1) in EV:
        d.setdefault(m.group(1), {})[int(m.group(2))] = (float(m.group(4)), float(m.group(5)))
for e in EV:
    if e not in d: continue
    r = d[e]; t = [r[k][0] for k in sorted(r)]; f = [r[k][1] for k in sorted(r)]
    imax = max(range(len(t)), key=lambda i: t[i]); fmax = max(range(len(f)), key=lambda i: f[i])
    tm = sorted(t)[len(t)//2]
    print(f"{e:20s} nranks={len(t)} tmax={max(t):.3g}(rank {imax}) tmin={min(t):.3g} tmed={tm:.3g} max/min={max(t)/max(min(t),1e-30):.2f} max/med={max(t)/tm:.2f} flopmax_rank={fmax} flop max/min={max(f)/max(min(f),1):.2f}")
    if '-v' in sys.argv[0:1] or len(t) <= 8: print("   t:", " ".join(f"{x:.3g}" for x in t)); print("   f:", " ".join(f"{x:.3g}" for x in f))
