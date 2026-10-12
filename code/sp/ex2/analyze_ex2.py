"""Parse PETSc ex2 '-efield_monitor' lines: E: t sum Enorm lgEnorm Emax lgEmax chargesum m0 m1 m2 ||E||_L2 Np (step)
Energy: W = 0.5*m2 (kinetic, sum w v^2) + 0.5*||E||^2 (field).  Landau rate: fit of log|E|_L2 peaks."""
import numpy as np, sys, json
def load(fn):
    rows = []
    for line in open(fn):
        if line.startswith("E: "):
            p = line.split()
            rows.append([float(p[1]), float(p[5]), float(p[7]), float(p[8]), float(p[9]), float(p[10]), float(p[11]), float(p[12])])
    a = np.array(rows)
    return dict(t=a[:, 0], Emax=a[:, 1], charge=a[:, 2], m0=a[:, 3], m1=a[:, 4], m2=a[:, 5], EL2=a[:, 6], Np=int(a[0, 7]))
def analyze(fn, tfit=(0.5, 15.5)):
    d = load(fn)
    W = 0.5 * d["m2"] + 0.5 * d["EL2"] ** 2
    Ek = 0.5 * d["m2"]
    t, e = d["t"], d["EL2"]
    # envelope peaks: max of |E| in consecutive windows of half the Langmuir period pi/omega_r (omega_r=1.4156 for k=0.5)
    win = np.pi / 1.4156; pk = []
    w0 = tfit[0]
    while w0 + win <= tfit[1]:
        idx = np.where((t >= w0) & (t < w0 + win))[0]
        if len(idx): pk.append(idx[np.argmax(e[idx])])
        w0 += win
    gamma = float(np.polyfit(t[pk], np.log(e[pk]), 1)[0]) if len(pk) >= 3 else float("nan")
    return dict(file=fn, Np=d["Np"], steps=len(t), t_end=float(t[-1]),
                W0=float(W[0]), W_rel_drift_max=float(np.max(np.abs(W - W[0])) / W[0]),
                W_rel_drift_end=float((W[-1] - W[0]) / W[0]),
                field_energy_max=float(np.max(0.5 * e ** 2)),
                charge_rel_drift_max=float(np.max(np.abs(d["charge"] - d["charge"][0])) / abs(d["charge"][0])),
                momentum_abs_max=float(np.max(np.abs(d["m1"]))),
                landau_gamma_fit=gamma, landau_gamma_theory=-0.1533, n_peaks=len(pk)), d, W
if __name__ == "__main__":
    for fn in sys.argv[1:]:
        r, _, _ = analyze(fn); print(json.dumps(r))
