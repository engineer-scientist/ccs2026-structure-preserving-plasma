"""2D ideal induction equation dB/dt = curl(u x B) on a periodic box.
CT scheme  : vector potential A_z on cell corners, B = curl A on faces (staggered, Yee/Evans-Hawley CT).
             Discrete div(curl) == 0 identically -> div B stays at round-off for all time.
Baseline   : cell-centred Bx, By, each advected with the same 2nd-order limited upwind scheme
             (non-CT, component-wise) -> div B is not controlled.
Velocity   : uniform diagonal (field-loop test, Gardiner & Stone 2005) or solid-body rotation.
Time       : SSP-RK2 (Heun).  Author: Grok Bot for S. Sharma, 2026-10-12.
"""
import numpy as np, json, sys

def minmod(a, b):
    return np.where(a * b > 0, np.sign(a) * np.minimum(np.abs(a), np.abs(b)), 0.0)

def adv_rhs(q, ux, uy, dx, dy, limiter=True):
    """-(u.grad) q, 2nd-order MUSCL upwind, periodic; ux,uy at the same points as q."""
    out = np.zeros_like(q)
    for axis, u, h in ((0, ux, dx), (1, uy, dy)):
        qm = np.roll(q, 1, axis); qp = np.roll(q, -1, axis)
        if limiter:
            s = minmod(q - qm, qp - q)
        else:
            s = 0.5 * (qp - qm)  # Fromm (unlimited, linear)
        # face i+1/2 states
        uf = 0.5 * (u + np.roll(u, -1, axis))
        qL = q + 0.5 * s
        qR = np.roll(q - 0.5 * s, -1, axis)
        flux_q = np.where(uf > 0, qL, qR)          # upwind face value
        # non-conservative advection: u * dq/dx ~ u * (q_{i+1/2} - q_{i-1/2}) / h
        out -= u * (flux_q - np.roll(flux_q, 1, axis)) / h
    return out

def setup(N, vel, case):
    L = 1.0 if vel == "uniform" else 2.0
    dx = dy = L / N
    xc = (np.arange(N) + 0.5) * dx - (0 if vel == "uniform" else 1.0)
    xn = np.arange(N) * dx - (0 if vel == "uniform" else 1.0)
    return L, dx, dy, xc, xn

def A_init(X, Y, case, vel):
    if case == "loop":
        x0, y0 = (0.5, 0.5) if vel == "uniform" else (0.4, 0.0)
        r = np.sqrt((X - x0) ** 2 + (Y - y0) ** 2); R0 = 0.3 if vel == "uniform" else 0.3
        return 1e-3 * np.maximum(R0 - r, 0.0)
    else:  # smooth periodic
        return np.sin(2 * np.pi * X) * np.sin(2 * np.pi * Y) / (2 * np.pi)

def vel_field(X, Y, vel):
    if vel == "uniform":
        return np.full_like(X, 1.0), np.full_like(Y, 0.5)
    return -Y * np.pi, X * np.pi  # solid body rotation, period 2

def B_from_A_faces(A, dx, dy):
    # A at corners (i,j) = (x_i, y_j). Bx at x-faces (x_i, y_{j+1/2}) = dA/dy ; By at y-faces (x_{i+1/2}, y_j) = -dA/dx
    Bx = (np.roll(A, -1, 1) - A) / dy
    By = -(np.roll(A, -1, 0) - A) / dx
    return Bx, By

def divB_faces(Bx, By, dx, dy):
    # cell (i+1/2, j+1/2): Bx faces at i and i+1 (index [i,j]), By faces at j and j+1 (index [i,j])
    return (np.roll(Bx, -1, 0) - Bx) / dx + (np.roll(By, -1, 1) - By) / dy

def faces_to_centres(Bx, By):
    return 0.5 * (Bx + np.roll(Bx, -1, 0)), 0.5 * (By + np.roll(By, -1, 1))

def divB_central(Bx, By, dx, dy):
    return (np.roll(Bx, -1, 0) - np.roll(Bx, 1, 0)) / (2 * dx) + (np.roll(By, -1, 1) - np.roll(By, 1, 1)) / (2 * dy)

def run(N=128, vel="uniform", case="loop", T=None, cfl=0.4, nrec=200, limiter=True):
    L, dx, dy, xc, xn = setup(N, vel, case)
    Xn, Yn = np.meshgrid(xn, xn, indexing="ij")          # corners
    Xc, Yc = np.meshgrid(xc, xc, indexing="ij")          # centres
    A = A_init(Xn, Yn, case, vel)
    uxn, uyn = vel_field(Xn, Yn, vel)
    uxc, uyc = vel_field(Xc, Yc, vel)
    umax = max(np.abs(uxn).max(), np.abs(uyn).max())
    dt = cfl * dx / umax
    if T is None:
        T = 2.0 if vel == "uniform" else 2.0          # uniform: 2 crossings in x ; rotation: one turn
    nsteps = int(np.ceil(T / dt)); dt = T / nsteps
    # baseline: cell-centred B obtained from the same initial A (div-free at t=0 in central sense ~ O(h^2))
    Bxf, Byf = B_from_A_faces(A, dx, dy)
    # baseline initial data: central-difference curl of A sampled at cell centres -> central div B == 0 at t=0
    Ac = A_init(Xc, Yc, case, vel)
    Bxc = (np.roll(Ac, -1, 1) - np.roll(Ac, 1, 1)) / (2 * dy)
    Byc = -(np.roll(Ac, -1, 0) - np.roll(Ac, 1, 0)) / (2 * dx)
    Bmax0 = max(np.abs(Bxf).max(), np.abs(Byf).max())
    hist = {"t": [], "ct_divB": [], "cc_divB": [], "ct_EB": [], "cc_EB": []}
    rec_every = max(1, nsteps // nrec)
    def record(t):
        bx, by = B_from_A_faces(A, dx, dy)
        hist["t"].append(t)
        hist["ct_divB"].append(float(np.abs(divB_faces(bx, by, dx, dy)).max() * dx / Bmax0))
        hist["cc_divB"].append(float(np.abs(divB_central(Bxc, Byc, dx, dy)).max() * dx / Bmax0))
        hist["ct_EB"].append(float(0.5 * np.sum(bx ** 2 + by ** 2) * dx * dy))
        hist["cc_EB"].append(float(0.5 * np.sum(Bxc ** 2 + Byc ** 2) * dx * dy))
    record(0.0)
    for n in range(1, nsteps + 1):
        # SSP-RK2 for A (CT) and for (Bxc, Byc) (baseline)
        A1 = A + dt * adv_rhs(A, uxn, uyn, dx, dy, limiter)
        A = 0.5 * (A + A1 + dt * adv_rhs(A1, uxn, uyn, dx, dy, limiter))
        if vel == "uniform":
            # uniform u: induction eq reduces to pure advection of each component
            bx1 = Bxc + dt * adv_rhs(Bxc, uxc, uyc, dx, dy, limiter)
            by1 = Byc + dt * adv_rhs(Byc, uxc, uyc, dx, dy, limiter)
            Bxc = 0.5 * (Bxc + bx1 + dt * adv_rhs(bx1, uxc, uyc, dx, dy, limiter))
            Byc = 0.5 * (Byc + by1 + dt * adv_rhs(by1, uxc, uyc, dx, dy, limiter))
        else:
            # rotation u=(-wy, wx): dB/dt = -(u.grad)B + (B.grad)u ; (B.grad)u = w(-By, Bx)
            w = np.pi
            def rhs(bx, by):
                return (adv_rhs(bx, uxc, uyc, dx, dy, limiter) - w * by,
                        adv_rhs(by, uxc, uyc, dx, dy, limiter) + w * bx)
            r1 = rhs(Bxc, Byc); bx1 = Bxc + dt * r1[0]; by1 = Byc + dt * r1[1]
            r2 = rhs(bx1, by1)
            Bxc = 0.5 * (Bxc + bx1 + dt * r2[0]); Byc = 0.5 * (Byc + by1 + dt * r2[1])
        if n % rec_every == 0 or n == nsteps:
            record(n * dt)
    # error vs exact (exact = initial for uniform after integer periods; rotation after full turn)
    A_ex = A_init(Xn, Yn, case, vel)
    bx, by = B_from_A_faces(A, dx, dy); bxe, bye = B_from_A_faces(A_ex, dx, dy)
    err_ct = float(np.sqrt(np.mean((bx - bxe) ** 2 + (by - bye) ** 2)))
    Ace = A_init(Xc, Yc, case, vel)
    bxce = (np.roll(Ace, -1, 1) - np.roll(Ace, 1, 1)) / (2 * dy); byce = -(np.roll(Ace, -1, 0) - np.roll(Ace, 1, 0)) / (2 * dx)
    err_cc = float(np.sqrt(np.mean((Bxc - bxce) ** 2 + (Byc - byce) ** 2)))
    return dict(N=N, vel=vel, case=case, T=T, dt=dt, nsteps=nsteps, err_ct=err_ct, err_cc=err_cc, hist=hist)

if __name__ == "__main__":
    out = sys.argv[1] if len(sys.argv) > 1 else "ct_results.json"
    res = {}
    res["loop_uniform"] = run(128, "uniform", "loop", T=2.0)
    res["loop_rotation"] = run(128, "rotation", "loop", T=2.0)
    conv = []
    for N in (32, 64, 128, 256):
        r = run(N, "uniform", "smooth", T=2.0, limiter=False, nrec=5)
        conv.append(dict(N=N, err_ct=r["err_ct"], err_cc=r["err_cc"], ct_divB_max=max(r["hist"]["ct_divB"])))
    res["convergence"] = conv
    json.dump(res, open(out, "w"))
    for k in ("loop_uniform", "loop_rotation"):
        h = res[k]["hist"]
        print(k, "cc divB t0", h["cc_divB"][0], "steps", res[k]["nsteps"], "max CT divB", max(h["ct_divB"]), "max CC divB", max(h["cc_divB"]),
              "EB ct", h["ct_EB"][0], "->", h["ct_EB"][-1], "cc", h["cc_EB"][0], "->", h["cc_EB"][-1])
    for c in conv: print(c)
