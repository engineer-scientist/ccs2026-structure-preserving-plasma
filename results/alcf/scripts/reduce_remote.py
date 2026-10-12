# Pure-python reducer (no numpy) for PETSc ex2 '-efield_monitor' output + -log_view; prints one JSON line per run.
import sys, glob, json, math, re, os
def load(fn):
    t=[];m1=[];m2=[];el2=[];np_=None
    for line in open(fn):
        if line.startswith("E: "):
            p=line.split(); t.append(float(p[1])); m1.append(float(p[9])); m2.append(float(p[10])); el2.append(float(p[11])); np_=int(p[12])
    return t,m1,m2,el2,np_
def logview(fn):
    out={}
    if not os.path.exists(fn): return out
    for line in open(fn):
        for ev in ("TSStep","DMSwarmMigrate","SNESSolve","KSPSolve"):
            if line.startswith(ev+" "):
                out[ev]=float(line.split()[3])
        if line.startswith("Time (sec):"): out["total"]=float(line.split()[2])
    return out
def fitgamma(t,e,t0=0.5,t1=15.5,win=math.pi/1.4156):
    xs=[];ys=[];w=t0
    while w+win<=t1:
        idx=[i for i in range(len(t)) if w<=t[i]<w+win]
        if idx:
            i=max(idx,key=lambda k:e[k]); xs.append(t[i]); ys.append(math.log(e[i]))
        w+=win
    n=len(xs)
    if n<3: return float('nan')
    mx=sum(xs)/n; my=sum(ys)/n
    return sum((x-mx)*(y-my) for x,y in zip(xs,ys))/sum((x-mx)**2 for x in xs)
for d in sys.argv[1:]:
    for fn in sorted(glob.glob(d+"/ex2_*.out")):
        t,m1,m2,el2,np_=load(fn)
        if not t: print(json.dumps({"file":fn,"error":"no E lines"})); continue
        W=[0.5*a+0.5*b*b for a,b in zip(m2,el2)]
        tag=os.path.basename(fn)[4:-4]
        r={"dir":os.path.basename(d),"run":tag,"Np":np_,"steps":len(t),"t_end":t[-1],
           "W_rel_drift_max":max(abs(w-W[0]) for w in W)/W[0],
           "W_drift_over_fieldE_max":max(abs(w-W[0]) for w in W)/max(0.5*b*b for b in el2),
           "momentum_abs_max":max(abs(x) for x in m1),
           "gamma":fitgamma(t,el2),"logview":logview(d+"/logview_"+tag+".txt"),
           "series_t":t[::10],"series_EL2":el2[::10],"series_W":[(w-W[0])/W[0] for w in W[::10]]}
        print(json.dumps(r))
