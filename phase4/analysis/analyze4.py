import glob, os, re, json, csv, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
D="/workspace/research/phase4/data"; F="/workspace/research/results/figs"; out={}
def summ(d):
    r={}
    for l in open(d+"/summary.txt"):
        p=l.split()
        if len(p)>2 and "=" in p[1]: r[p[0]]=dict(x.split("=",1) for x in p[1:] if "=" in x)
    return r
def tsstep(f):
    for l in open(f):
        if l.startswith("TSStep "): return float(l.split()[3])
# strong
fig,ax=plt.subplots(1,2,figsize=(11,4.2)); rows=[]
for m,d,ppn,c,mk in (("Aurora",glob.glob(D+"/p4_results_aurora/runs_*ss_*")[0],96,"C1","s"),("Polaris",glob.glob(D+"/p4_results_polaris/runs_*ss_*")[0],32,"C0","o")):
    s=summ(d); n=sorted(int(v["np"]) for v in s.values()); nodes=np.array(n)//ppn
    w=np.array([float(s[f"ss_np{k}"]["wall_s"]) for k in n]); t=np.array([tsstep(f"{d}/logview_ss_np{k}.txt") for k in n])
    for k,nd,a,b in zip(n,nodes,w,t): rows.append(dict(machine=m,ranks=k,nodes=int(nd),wall_s=a,TSStep_s=b,speedup_wall=w[0]/a,eff_wall=w[0]/a/nd,speedup_ts=t[0]/b,eff_ts=t[0]/b/nd,Np_conserved=s[f"ss_np{k}"]["Np_first"]==s[f"ss_np{k}"]["Np_last"]))
    ax[0].loglog(nodes,w,mk+"-",color=c,label=f"{m} wall"); ax[0].loglog(nodes,w[0]/nodes,":",color=c,alpha=.6)
    ax[1].semilogx(nodes,w[0]/w/nodes,mk+"-",color=c,label=f"{m} (wall)"); ax[1].semilogx(nodes,t[0]/t/nodes,mk+"--",color=c,mfc="none",label=f"{m} (TSStep only)")
ax[0].set_xlabel("nodes"); ax[0].set_ylabel("wall time, 100 steps (s)"); ax[0].set_title("Strong scaling: fixed 9.83M particles, 49,152 cells (dotted = ideal)",fontsize=9); ax[0].legend(fontsize=8)
ax[1].axhline(1,color="gray",ls=":"); ax[1].set_ylim(0,1.1); ax[1].set_xlabel("nodes"); ax[1].set_ylabel("strong-scaling efficiency vs 1 node"); ax[1].legend(fontsize=7); ax[1].set_title("CPU path, 96 (Aurora) / 32 (Polaris) ranks per node",fontsize=9)
plt.tight_layout(); plt.savefig(F+"/strong_scaling.png",dpi=150); plt.close()
with open(D+"/strong_scaling.csv","w") as f: wr=csv.DictWriter(f,rows[0].keys()); wr.writeheader(); wr.writerows(rows)
out["strong"]=rows
# weak
p3=list(csv.DictReader(open("/workspace/research/phase3/data/weak_scaling_v2.csv")))
fig,ax=plt.subplots(figsize=(7,4.3))
for var,mach,lab,c,mk in (("after_fix","aurora","Aurora phase 3 (1-16 nodes)","C1","s"),("after_fix","polaris","Polaris phase 3 (1-8 nodes)","C0","o"),("before_growbox","aurora","Aurora before fix","C1","x"),("before_growbox","polaris","Polaris before fix","C0","x")):
    r=[x for x in p3 if x["variant"]==var and x["machine"]==mach and x["eff_vs_1node"]]
    ax.semilogx([int(x["nodes"]) for x in r],[float(x["eff_vs_1node"]) for x in r],mk+("-" if var=="after_fix" else ":"),color=c,label=lab,alpha=1 if var=="after_fix" else .5)
wsd={}
for d in glob.glob(D+"/p4_results_aurora/runs_*ws*"):
    s=summ(d)
    for k,v in s.items(): wsd[int(v["np"])]=(float(v["wall_s"]),tsstep(f"{d}/logview_{k}.txt"),v["Np_first"]==v["Np_last"],int(v["particles"]))
w96,t96=wsd[96][0],wsd[96][1]; ws=[]
for k in sorted(wsd): ws.append(dict(ranks=k,nodes=k//96,wall_s=wsd[k][0],TSStep_s=wsd[k][1],eff_wall=w96/wsd[k][0],eff_ts=t96/wsd[k][1],particles=wsd[k][3],Np_conserved=wsd[k][2]))
ax.semilogx([x["nodes"] for x in ws if x["nodes"]>1],[x["eff_ts"] for x in ws if x["nodes"]>1],"D-",color="C3",ms=8,label="Aurora phase 4: 32 & 64 nodes")
ax.axhline(1,color="gray",ls=":"); ax.set_ylim(0,1.1); ax.set_xlabel("nodes (51,200 particles per rank)"); ax.set_ylabel("weak-scaling efficiency vs 1 node (TSStep)"); ax.legend(fontsize=7); ax.set_title("ex2 weak scaling with the wrap-cell hash fix, up to 6,144 ranks / 315M particles",fontsize=9)
plt.tight_layout(); plt.savefig(F+"/weak_scaling_v3.png",dpi=150); plt.close()
with open(D+"/weak_scaling_p4.csv","w") as f: wr=csv.DictWriter(f,ws[0].keys()); wr.writeheader(); wr.writerows(ws)
out["weak"]=ws; out["p3_aurora_np96_TSStep_s"]=float([x for x in p3 if x["variant"]=="after_fix" and x["machine"]=="aurora" and x["ranks"]=="96"][0]["s_per_step"])*100
# landau
fig,ax=plt.subplots(2,3,figsize=(14,7)); L=[]
for j,run in enumerate(["x10_dt05","x10_dt01","x100_dt05_lowres"]):
  ser={}
  for m,c in (("aurora","C1"),("polaris","C0")):
    f=glob.glob(f"{D}/p4_results_{m}/runs_*lr_*/ex2_{run}.out")[0]
    E=np.array([[float(e.split()[1]),float(e.split()[10]),float(e.split()[11]),float(e.split()[12])] for e in open(f,errors="replace") if e.startswith("E: ")])
    t=E[:,0]; W=0.5*E[:,1]+0.5*E[:,2]**2; rel=(W-W[0])/abs(W[0]); EL=E[:,2]; ser[m]=(t,W,EL)
    n=len(t); q=n//4
    # damping fit on local maxima in t in [0.5,15.5]
    pk=[i for i in range(1,n-1) if EL[i]>=EL[i-1] and EL[i]>=EL[i+1] and 0.5<=t[i]<=15.5]
    g=np.polyfit(t[pk],np.log(EL[pk]),1)[0] if len(pk)>2 else float("nan")
    # secular test: linear fit of rel over second half
    h=n//2; slope=np.polyfit(t[h:],rel[h:],1)[0]
    L.append(dict(machine=m,run=run,steps=n,t_end=t[-1],Np0=int(E[0,3]),Np1=int(E[-1,3]),max_abs_rel=float(np.max(np.abs(rel))),rel_end=float(rel[-1]),mean_rel_q1=float(rel[:q].mean()),mean_rel_q4=float(rel[-q:].mean()),max_q1=float(np.abs(rel[:q]).max()),max_q4=float(np.abs(rel[-q:]).max()),slope_2nd_half_per_t=float(slope),gamma_fit=float(g),EL2_min=float(EL.min()),EL2_end=float(EL[-1])))
    ax[0,j].semilogy(t,EL,color=c,lw=.6,label=m)
    ax[1,j].plot(t,rel,color=c,lw=.6,label=m)
  tt=np.linspace(0,min(40,ser["aurora"][0][-1]),50); ax[0,j].semilogy(tt,ser["aurora"][2][0]*1.5*np.exp(-0.1533*tt),"k--",lw=.8,label="theory γ=-0.1533")
  a,p=ser["aurora"],ser["polaris"]; L[-1]["aurora_vs_polaris_max_dW_rel"]=float(np.max(np.abs(a[1]-p[1]))/abs(a[1][0])); L[-1]["aurora_vs_polaris_max_dEL2_rel"]=float(np.max(np.abs(a[2]-p[2])/np.abs(p[2])))
  ax[0,j].set_title(f"{run}: t_end={a[0][-1]:.0f}, {len(a[0])} steps",fontsize=9); ax[0,j].set_ylabel(r"$\|E\|_{L^2}$"); ax[0,j].legend(fontsize=7)
  ax[1,j].set_xlabel("t"); ax[1,j].set_ylabel(r"$(W-W_0)/W_0$"); ax[1,j].legend(fontsize=7)
plt.suptitle("Long-run Landau damping (PETSc ex2, TSBASICSYMPLECTIC), Aurora vs Polaris",fontsize=10); plt.tight_layout(); plt.savefig(F+"/landau_longrun.png",dpi=150); plt.close()
out["landau"]=L
# pushers
B="/workspace/research/phase4/box_pusher"; fig,ax=plt.subplots(2,2,figsize=(11,7)); P=[]
for i,fld in enumerate(["uniform","solovev"]):
  for m,c in (("boris","C2"),("rk4","C3"),("dg","C4")):
    a=np.array([list(map(float,l.split())) for l in open(f"{B}/{fld}_{m}_long.dat") if len(l.split())==6 and l.endswith("\n")]); t=a[:,0]; s=(a[:,0]/0.2).astype(np.int64)
    P.append(dict(field=fld,method=m,last_step=int(round(t[-1]/0.2)),complete=bool(round(t[-1]/0.2)>=1e8),max_abs_dH=float(np.abs(a[:,1]).max()),dH_end=float(a[-1,1]),max_abs_dP=float(np.abs(a[:,2]).max()),dP_end=float(a[-1,2])))
    lab=m+("" if P[-1]["complete"] else f" (stalled: SNES failures cut Δt, reached t={t[-1]:.0f} of 2e7)")
    ax[0,i].loglog(np.maximum(s[1:],1),np.abs(a[1:,1])+1e-17,color=c,lw=.7,label=lab); ax[1,i].loglog(np.maximum(s[1:],1),np.abs(a[1:,2])+1e-17,color=c,lw=.7,label=lab)
  ax[0,i].set_title(f"{fld} B: |ΔH/H0|",fontsize=9); ax[1,i].set_title(f"{fld} B: |Δp/p0| ("+("p_par" if fld=="uniform" else "canonical p_φ")+")",fontsize=9)
  for k in (0,1): ax[k,i].set_xlabel("step (Δt=0.2)"); ax[k,i].legend(fontsize=7)
plt.suptitle("Charged-particle pushers, 10^8 steps (box)",fontsize=10); plt.tight_layout(); plt.savefig(F+"/pusher_longrun.png",dpi=150); plt.close()
out["pusher"]=P
json.dump(out,open(D+"/phase4_analysis.json","w"),indent=1,default=str); print(json.dumps(out,indent=1,default=str))
