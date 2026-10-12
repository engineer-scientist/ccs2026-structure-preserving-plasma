import numpy as np, matplotlib
def ld(fn):
    rows=[l.split() for l in open(fn)]
    return np.array([[float(x) for x in r] for r in rows if len(r)==6], ndmin=2)
import matplotlib, json, os; matplotlib.use("Agg"); import matplotlib.pyplot as plt
R = "/workspace/research/results/pusher"; F = "/workspace/research/results/figs"
lab = {"boris": "Boris (volume-preserving)", "rk4": "TSRK 4 (RK4)", "bsi": "TSBASICSYMPLECTIC (Verlet)", "dg": "TSDISCGRAD (Gonzalez)"}
col = {"boris": "C0", "rk4": "C3", "bsi": "C2", "dg": "C1"}
summary = {}
for field, title, dt in (("es", "Electrostatic (B=0, pendulum+harmonic)", 0.1), ("uniform", r"Uniform $B_z$ + radial E well ($\omega_c\Delta t=0.2$)", 0.2), ("solovev", r"Solov'ev tokamak field ($\omega_c\Delta t\approx0.2$)", 0.2)):
    fig, ax = plt.subplots(1, 2 if field != "es" else 1, figsize=(11 if field != "es" else 6, 4)); ax = np.atleast_1d(ax)
    for m in ("rk4", "bsi", "boris", "dg"):
        fn = f"{R}/{field}_{m}.dat"
        if not os.path.exists(fn) or os.path.getsize(fn) == 0: continue
        a = ld(fn); t = a[:, 0]; dH = np.abs(a[:, 1]); dP = np.abs(a[:, 2])
        ok = np.isfinite(dH) & (dH < 1e3)
        summary[f"{field}_{m}"] = dict(t_end=float(t[-1]), steps=int(round(t[-1] / dt)), energy_rel_err_max=float(np.nanmax(np.where(np.isfinite(dH), dH, np.inf))), energy_rel_err_end=float(dH[-1]), pphi_rel_err_max=float(np.nanmax(np.where(np.isfinite(dP), dP, np.inf))), blew_up=bool(not ok.all()))
        tt = t[ok][1:]; ax[0].loglog(tt, np.maximum(dH[ok][1:], 1e-17), color=col[m], label=lab[m] + (" (blows up)" if not ok.all() else ""))
        if field != "es": ax[1].loglog(tt, np.maximum(dP[ok][1:], 1e-17), color=col[m], label=lab[m])
    ax[0].set_xlabel("time"); ax[0].set_ylabel(r"$|H(t)-H_0|/|H_0|$"); ax[0].set_title("Energy error: " + title, fontsize=9); ax[0].legend(fontsize=7); ax[0].set_ylim(1e-17, 10)
    if field != "es":
        ax[1].set_xlabel("time"); ax[1].set_ylabel(r"$|p_\varphi(t)-p_{\varphi,0}|/|p_{\varphi,0}|$"); ax[1].set_title(r"Canonical angular momentum $p_\varphi$ error", fontsize=9); ax[1].legend(fontsize=7); ax[1].set_ylim(1e-17, 10)
    plt.tight_layout(); plt.savefig(f"{F}/pusher_{field}.png", dpi=150); plt.close()
# orbit picture (Solov'ev): R-Z projection
fig, ax = plt.subplots(figsize=(4, 5))
for m in ("rk4", "boris", "dg"):
    fn = f"{R}/solovev_{m}.dat"
    if os.path.exists(fn):
        a = ld(fn); Rr = np.hypot(a[:, 3], a[:, 4]); ax.plot(Rr, a[:, 5], ".", ms=1, color=col[m], label=lab[m])
ax.set_xlabel("R"); ax.set_ylabel("Z"); ax.set_title("Guiding orbit (R,Z), Solov'ev", fontsize=9); ax.legend(fontsize=6, markerscale=8)
plt.tight_layout(); plt.savefig(f"{F}/pusher_solovev_orbit.png", dpi=150)
json.dump(summary, open(f"{R}/summary.json", "w"), indent=1)
for k, v in summary.items(): print(k, v)
