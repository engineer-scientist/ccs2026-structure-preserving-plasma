# Phase-3 reducer: one CSV row per run dir/run. usage: python3 reduce3.py runs_dir [...]
import sys, glob, os, re
EVS = ["TSStep","DMSwarmMigrate","DMPlexLocatePoints","SNESSolve","DMCreateMassMat","RhsV"]
print("dir,run,np,ppn,nx,fix,steps_req,rc,Np0,Np1,Q0,Q1,W_rel_drift_max,steps_done," + ",".join(f"{e}_max,{e}_ratio" for e in EVS) + ",total_s,wall_s")
for d in sys.argv[1:]:
    summ = {}
    if os.path.exists(d+"/summary.txt"):
        for l in open(d+"/summary.txt"):
            p = l.split()
            if len(p) > 3 and "=" in p[1]:
                summ[p[0]] = dict(x.split("=",1) for x in p[1:] if "=" in x)
    for fn in sorted(glob.glob(d+"/ex2_*.out")):
        tag = os.path.basename(fn)[4:-4]; s = summ.get(tag, {})
        E = [l.split() for l in open(fn, errors="replace") if l.startswith("E: ")]
        Np0 = E[0][12] if E else ""; Np1 = E[-1][12] if E else ""
        Q0 = E[0][8] if E else ""; Q1 = E[-1][8] if E else ""  # global total weight (moment 0); col 7 is rank-0-local only
        W = [0.5*float(e[10]) + 0.5*float(e[11])**2 for e in E]
        drift = max(abs(w-W[0]) for w in W)/abs(W[0]) if W else ""
        lv = {}
        lf = d+"/logview_"+tag+".txt"
        tot = ""
        if os.path.exists(lf):
            for l in open(lf):
                p = l.split()
                if p and p[0] in EVS and len(p) > 4 and p[0] not in lv: lv[p[0]] = (p[3], p[4])
                if l.startswith("Time (sec):"): tot = p[2]
        row = [os.path.basename(d), tag, s.get("np",""), s.get("ppn",""), s.get("nx",""), s.get("fix",""), s.get("steps",""), s.get("rc",""), Np0, Np1, Q0, Q1, f"{drift:.4e}" if drift != "" else "", str(len(E))]
        for e in EVS: row += list(lv.get(e, ("","")))
        row += [tot, s.get("wall_s","")]
        print(",".join(row))
