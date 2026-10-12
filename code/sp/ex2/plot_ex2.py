import json, glob, os, numpy as np, matplotlib; matplotlib.use("Agg"); import matplotlib.pyplot as plt
A = "/workspace/research/results/alcf"; F = "/workspace/research/results/figs"
P = json.load(open(f"{A}/polaris_1node_7737974_summary.json")); Au = json.load(open(f"{A}/aurora_1node_8921878_summary.json"))
box = [json.loads(l) for l in open(f"{A}/box_landau.jsonl")][0]
fig, ax = plt.subplots(1, 2, figsize=(11, 4))
t = np.array(P["np1_series_t"]); ax[0].semilogy(t, P["np1_series_EL2"], "o-", ms=3, label="Polaris & Aurora, 1 rank, 51,200 particles (identical)")
ax[0].semilogy(box["series_t"], box["series_EL2"], "-", lw=0.8, alpha=0.7, label="box, 4,000 particles")
tt = np.linspace(0, 30, 50); ax[0].semilogy(tt, 3.0e-2 * np.exp(-0.1533 * tt), "k--", label=r"linear theory $\gamma=-0.1533$")
ax[0].set_xlabel("t"); ax[0].set_ylabel(r"$\|E\|_{L^2}$"); ax[0].set_title("Landau damping (k=0.5, α=0.01): PETSc ex2, DMSwarm + TSBASICSYMPLECTIC", fontsize=8); ax[0].legend(fontsize=7)
ax[1].plot(t, np.array(P["np1_series_Wrel"]) * 1e5, "o-", ms=3, label="1 rank")
ax[1].set_xlabel("t"); ax[1].set_ylabel(r"$(W-W_0)/W_0\ \times 10^{5}$"); ax[1].set_title(r"Total energy $W=\frac{1}{2}\sum_p w_p v_p^2+\frac{1}{2}\|E\|^2$ drift (same for 1-16 ranks)", fontsize=8)
plt.tight_layout(); plt.savefig(f"{F}/ex2_landau_energy.png", dpi=150); plt.close()
# weak scaling (pilot 4: dt=0.01, 100 steps, no particle loss)
fig, ax = plt.subplots(1, 2, figsize=(11, 4.2))
for name, files, mk, c in (("Polaris (CPU path, 32 ranks/node)", ["polaris_p4_2node_7738042"], "o", "C0"), ("Aurora (CPU path, 96 ranks/node)", ["aurora_p4_1node_8921920", "aurora_p4_2node_8921924"], "s", "C1")):
    rows = []
    for f in files:
        fn = f"{A}/{f}_summary.json"
        if not os.path.exists(fn): continue
        for r in json.load(open(fn))["runs"]:
            if "TSStep_s" in r and not (r["np"] == 32 and r.get("ppn") == 16) and not (r["np"] == 96 and r.get("ppn") == 48):
                rows.append((r["np"], r["TSStep_s"] / r["steps"], r["DMSwarmMigrate_s"] / r["steps"], r["SNESSolve_s"] / r["steps"], r.get("nodes", 1)))
    rows.sort(); n = np.array([r[0] for r in rows]); ts = np.array([r[1] for r in rows]); sn = np.array([r[3] for r in rows]); nd = np.array([r[4] for r in rows])
    ax[0].loglog(n, ts, mk + "-", color=c, label=name + ": TSStep")
    ax[0].loglog(n, sn, mk + ":", color=c, mfc="none", label=name + ": Poisson (SNESSolve)")
    ax[1].semilogx(n, ts[0] / ts, mk + "-", color=c, label=name)
    for a, b, k in zip(n, ts, nd):
        if k == 2: ax[1].annotate("2 nodes", (a, ts[0] / b), fontsize=7, xytext=(-10, 8), textcoords="offset points")
ax[0].set_xlabel("MPI ranks (51,200 particles + 256 cells per rank)"); ax[0].set_ylabel("seconds per time step"); ax[0].legend(fontsize=6)
ax[0].set_title("ex2 weak-scaling pilot (debug queue, <= 2 nodes)", fontsize=9)
ax[1].axhline(1, color="gray", ls=":"); ax[1].set_ylim(0, 1.1); ax[1].set_xlabel("MPI ranks"); ax[1].set_ylabel("weak-scaling efficiency  T(1)/T(P)")
ax[1].set_title("Particle location/migration (DMPlexLocatePoints) dominates; Poisson 1-2 %", fontsize=9); ax[1].legend(fontsize=7)
plt.tight_layout(); plt.savefig(f"{F}/ex2_weak_scaling_pilot.png", dpi=150)
print("ok")
