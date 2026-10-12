import csv, glob, os, sys
import matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D="/workspace/research/phase3/data"; R="/workspace/research/results"
def rows(f): return list(csv.DictReader(open(f)))
def ps(r): return float(r["TSStep_max"])/int(r["steps_done"] or r["steps_req"])
out=[]  # combined CSV rows
# Phase 2 (before, fixed box L=4pi, dt=0.01)
p2=rows(R+"/ex2_weak_scaling.csv")
for m in ("aurora","polaris"):
    rr=[r for r in p2 if r["machine"]==m and r["ranks_per_node"] in (r["ranks"],"96","32") ]
    rr=[r for r in rr if not (r["ranks"]=="96" and r["ranks_per_node"]=="48") and not (r["ranks"]=="32" and r["ranks_per_node"]=="16")]
    t1=float([r for r in rr if r["ranks"]=="1"][0]["TSStep_s_per_step"])
    for r in rr: out.append(dict(series=f"{m} phase2 before (fixed box, dt=0.01)",machine=m,variant="before_fixedbox",ranks=int(r["ranks"]),nodes=int(r["nodes"]),s_per_step=float(r["TSStep_s_per_step"]),eff_vs_1rank=t1/float(r["TSStep_s_per_step"]),locate_maxmin=r["LocatePoints_max_over_min"],Np_conserved="yes",W_rel_drift=r["W_rel_drift_max"]))
def add(f, m, tags):
    rr=rows(f)
    for tag,variant,label,ref in tags:
        sel=[r for r in rr if r["run"].startswith(tag) and r["rc"]=="0"]
        t1=ps([r for r in rr if r["run"]==ref][0]) if ref else None
        for r in sel:
            P=int(r["np"]); ppn=int(r["ppn"])
            out.append(dict(series=label,machine=m,variant=variant,ranks=P,nodes=max(1,P//ppn) if P>1 else 1,s_per_step=ps(r),eff_vs_1rank=(t1/ps(r)) if t1 else "",locate_maxmin=r["DMPlexLocatePoints_ratio"],Np_conserved="yes" if r["Np0"]==r["Np1"] and r["Q0"]==r["Q1"] else "NO",W_rel_drift=r["W_rel_drift_max"]))
A=D+"/aurora_p3ws_8922217_reduced.csv"
add(A,"aurora",[("ws_fix_np","after_fix","aurora phase3 after fix (grow box, dt=0.05)","ws_fix_np1"),("ws_base_np","before_growbox","aurora phase3 before fix (grow box, dt=0.05)",None)])
# base reference: unfixed np1 from diag (D6)
d6=[r for r in rows(D+"/aurora_p3diag_8922161_reduced.csv") if r["run"]=="D6_grow_np1_base"][0]
for o in out:
    if o["series"].startswith("aurora phase3 before"): o["eff_vs_1rank"]=ps(d6)/o["s_per_step"]
for f in glob.glob(D+"/polaris_p3ws_*_reduced.csv"):
    add(f,"polaris",[("ws_fix_np","after_fix","polaris phase3 after fix, CPU 32/node","ws_fix_np1"),("ws_base_np","before_growbox","polaris phase3 before fix, CPU 32/node","ws_fix_np1"),("gpu_fix_np","after_fix_gpu","polaris phase3 after fix, 4 ranks/node + CUDA vec/mat","ws_fix_np1"),("cpu4_fix_np","after_fix_cpu4","polaris phase3 after fix, 4 ranks/node CPU","ws_fix_np1")])
out.sort(key=lambda o:(o["series"],o["ranks"]))
for o in out:
    o["eff_vs_1node"]=""
for s in set(o["series"] for o in out):
    ss=[o for o in out if o["series"]==s and o["ranks"]>1]
    if ss:
        b=min(ss,key=lambda o:o["ranks"])
        for o in ss: o["eff_vs_1node"]=b["s_per_step"]/o["s_per_step"]
keys=["series","machine","variant","ranks","nodes","s_per_step","eff_vs_1rank","eff_vs_1node","locate_maxmin","Np_conserved","W_rel_drift"]
with open(D+"/weak_scaling_v2.csv","w",newline="") as fh:
    w=csv.DictWriter(fh,keys); w.writeheader(); [w.writerow({k:o[k] for k in keys}) for o in out]
# figure
fig,ax=plt.subplots(1,2,figsize=(13,5))
sty={"before_fixedbox":dict(ls=":",marker="x"),"before_growbox":dict(ls="--",marker="v"),"after_fix":dict(ls="-",marker="o",lw=2),"after_fix_gpu":dict(ls="-",marker="s"),"after_fix_cpu4":dict(ls="-.",marker="^")}
col={"aurora":"tab:blue","polaris":"tab:red"}
for s in sorted(set(o["series"] for o in out)):
    ss=sorted([o for o in out if o["series"]==s and o["eff_vs_1rank"]!=""],key=lambda o:o["ranks"])
    if not ss: continue
    v=ss[0]["variant"]; m=ss[0]["machine"]
    c=col[m] if v!="after_fix_gpu" else "tab:green"
    ax[0].semilogx([o["ranks"] for o in ss],[o["eff_vs_1rank"] for o in ss],color=c,label=s,**sty[v])
ax[0].set_xlabel("MPI ranks (51,200 particles + 256 cells per rank)"); ax[0].set_ylabel("weak-scaling efficiency  T(1 rank)/T(P)")
ax[0].set_ylim(0,1.1); ax[0].grid(alpha=.3); ax[0].legend(fontsize=7,loc="lower left"); ax[0].set_title("(a) Weak scaling, PETSc ex2 Vlasov-Poisson PIC, 1 to 16 nodes")
# breakdown aurora before vs after (per step)
ra=rows(A); Ps=[96,192,384,768,1536]; import numpy as np
x=np.arange(len(Ps)); wdt=0.2
def get(tag,P,k):
    r=[r for r in ra if r["run"]==f"{tag}{P}"][0]; return float(r[k])/int(r["steps_done"])
for i,(tag,lab,h) in enumerate([("ws_base_np","before",""),("ws_fix_np","after","//")]):
    ts=[get(tag,P,"TSStep_max") for P in Ps]; lo=[get(tag,P,"DMPlexLocatePoints_max") for P in Ps]; sn=[get(tag,P,"SNESSolve_max") for P in Ps]
    ax[1].bar(x+(2*i-1.5)*wdt*1.0,ts,wdt,color="0.5",hatch=h,label=f"TSStep (whole step), {lab} fix")
    ax[1].bar(x+(2*i-0.5)*wdt*1.0,lo,wdt,color="tab:orange",hatch=h,label=f"DMPlexLocatePoints (slowest rank), {lab} fix")
    ax[1].plot(x+(2*i-1)*wdt,sn,"k"+("s" if i else "v"),label=f"Poisson SNESSolve, {lab} fix")
ax[1].set_yscale("log"); ax[1].set_xticks(x); ax[1].set_xticklabels([f"{P}\n({P//96} node{'s' if P>96 else ''})" for P in Ps])
ax[1].set_ylabel("seconds per time step (max over ranks)"); ax[1].set_title("(b) Aurora time breakdown, before vs after wrap-rank hash fix")
ax[1].set_ylim(1e-2,300); ax[1].legend(fontsize=6.5,ncol=2,loc="upper center"); ax[1].grid(alpha=.3,axis="y")
fig.text(0.5,0.005,"Locate time includes the post-step migration outside TSStep, so it can exceed TSStep. Grow-the-box runs: L = pi*P, dt = 0.05, 100 steps (20 for 'before').",ha="center",fontsize=7)
plt.tight_layout(rect=(0,0.03,1,1)); plt.savefig("/workspace/research/phase3/figs/weak_scaling_v2.png",dpi=150)
for o in out: print(o["series"][:45].ljust(46),o["ranks"],f'{o["s_per_step"]:.4f}',f'{o["eff_vs_1node"]:.3f}' if o["eff_vs_1node"]!="" else "", f'{o["eff_vs_1rank"]:.3f}' if o["eff_vs_1rank"]!="" else "",o["locate_maxmin"],o["Np_conserved"])
