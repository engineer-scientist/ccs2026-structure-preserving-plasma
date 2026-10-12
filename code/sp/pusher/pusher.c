static char help[] = "Charged-particle pusher study: Boris vs PETSc TSRK4 / TSBASICSYMPLECTIC / TSDISCGRAD.\n"
                     "  -field es|uniform|solovev  -method boris|rk4|bsi|dg  -dt  -steps  -nout  -o file\n";
/* State u = (x,y,z,vx,vy,vz), normalized units, a = qm (E + v x B), E = -grad(phi).
   Invariants: H = |v|^2/2 + qm*phi(x); for axisymmetric 'solovev' field also
   p_phi = (x vy - y vx) + qm*psi(R,Z)  (canonical toroidal momentum, A_phi = psi/R). */
#include <petscts.h>

typedef enum { F_ES, F_UNIFORM, F_SOLOVEV } FieldType;
static const char *FieldTypes[] = {"es", "uniform", "solovev", "FieldType", "F_", NULL};
typedef struct {
  FieldType field;
  PetscReal qm, B0, kE, q0, kappa;
  PetscInt  nout;
  FILE     *fp;
  PetscReal H0, P0, maxdH, maxdP;
  Vec       U; /* full state, for BSI momentum substep (needs v for v x B) */
} Ctx;

static void Fields(Ctx *c, const PetscReal *x, PetscReal E[3], PetscReal B[3], PetscReal *phi, PetscReal *psi)
{
  E[0] = E[1] = E[2] = B[0] = B[1] = B[2] = 0;
  *phi = 0; *psi = 0;
  if (c->field == F_ES) { /* pendulum in x, harmonic in y,z, B=0 : separable Hamiltonian */
    *phi = (1.0 - PetscCosReal(x[0])) + 0.5 * (x[1] * x[1] + x[2] * x[2]);
    E[0] = -PetscSinReal(x[0]); E[1] = -x[1]; E[2] = -x[2];
  } else if (c->field == F_UNIFORM) { /* uniform Bz + radial harmonic electrostatic well */
    B[2] = c->B0;
    *phi = 0.5 * c->kE * (x[0] * x[0] + x[1] * x[1]);
    E[0] = -c->kE * x[0]; E[1] = -c->kE * x[1];
    *psi = 0.5 * c->B0 * (x[0] * x[0] + x[1] * x[1]); /* R*A_phi, A_phi = B0 R/2 */
  } else { /* Solov'ev equilibrium, R0 = B0 = 1 */
    PetscReal R = PetscSqrtReal(x[0] * x[0] + x[1] * x[1]), Z = x[2], cp = x[0] / R, sp = x[1] / R;
    PetscReal a = 1.0 / (2.0 * c->q0 * c->kappa), k2 = c->kappa * c->kappa;
    PetscReal BR = -2.0 * a * R * Z, BZ = a * (2.0 * Z * Z + k2 * (R * R - 1.0)), Bp = 1.0 / R;
    *psi = a * (R * R * Z * Z + 0.25 * k2 * (R * R - 1.0) * (R * R - 1.0));
    B[0] = BR * cp - Bp * sp; B[1] = BR * sp + Bp * cp; B[2] = BZ;
  }
}

static void Invariants(Ctx *c, const PetscReal *u, PetscReal *H, PetscReal *P)
{
  PetscReal E[3], B[3], phi, psi;
  Fields(c, u, E, B, &phi, &psi);
  *H = 0.5 * (u[3] * u[3] + u[4] * u[4] + u[5] * u[5]) + c->qm * phi;
  *P = (u[0] * u[4] - u[1] * u[3]) + c->qm * psi;
}

static void Record(Ctx *c, PetscInt step, PetscReal t, const PetscReal *u)
{
  if (!c->fp) return;
  PetscReal H, P;
  Invariants(c, u, &H, &P);
  if (step == 0) { c->H0 = H; c->P0 = P; c->maxdH = 0; c->maxdP = 0; }
  c->maxdH = PetscMax(c->maxdH, PetscAbsReal((H - c->H0) / c->H0));
  c->maxdP = PetscMax(c->maxdP, PetscAbsReal((P - c->P0) / (PetscAbsReal(c->P0) > 0 ? PetscAbsReal(c->P0) : 1.0)));
  if (step % c->nout == 0) fprintf(c->fp, "%.10e %.16e %.16e %.10e %.10e %.10e\n", (double)t, (double)((H - c->H0) / PetscAbsReal(c->H0)), (double)((P - c->P0) / (PetscAbsReal(c->P0) > 0 ? PetscAbsReal(c->P0) : 1.0)), (double)u[0], (double)u[1], (double)u[2]);
}

static PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal t, Vec U, void *ctx)
{
  const PetscScalar *u;
  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u));
  Record((Ctx *)ctx, step, t, u);
  PetscCall(VecRestoreArrayRead(U, &u));
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode RHS(TS ts, PetscReal t, Vec U, Vec F, void *ctx)
{
  Ctx *c = (Ctx *)ctx; const PetscScalar *u; PetscScalar *f; PetscReal E[3], B[3], phi, psi;
  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u)); PetscCall(VecGetArray(F, &f));
  Fields(c, u, E, B, &phi, &psi);
  f[0] = u[3]; f[1] = u[4]; f[2] = u[5];
  f[3] = c->qm * (E[0] + u[4] * B[2] - u[5] * B[1]);
  f[4] = c->qm * (E[1] + u[5] * B[0] - u[3] * B[2]);
  f[5] = c->qm * (E[2] + u[3] * B[1] - u[4] * B[0]);
  PetscCall(VecRestoreArrayRead(U, &u)); PetscCall(VecRestoreArray(F, &f));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* BSI split: position update x' = v ; momentum update v' = qm(E(x) + v x B(x)) with v taken from current state */
static PetscErrorCode RHSPos(TS ts, PetscReal t, Vec V, Vec Xres, void *ctx)
{
  PetscFunctionBeginUser;
  PetscCall(VecCopy(V, Xres));
  PetscFunctionReturn(PETSC_SUCCESS);
}
static PetscErrorCode RHSMom(TS ts, PetscReal t, Vec X, Vec Vres, void *ctx)
{
  Ctx *c = (Ctx *)ctx; const PetscScalar *x, *u; PetscScalar *f; PetscReal E[3], B[3], phi, psi;
  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(X, &x)); PetscCall(VecGetArrayRead(c->U, &u)); PetscCall(VecGetArray(Vres, &f));
  Fields(c, x, E, B, &phi, &psi);
  f[0] = c->qm * (E[0] + u[4] * B[2] - u[5] * B[1]);
  f[1] = c->qm * (E[1] + u[5] * B[0] - u[3] * B[2]);
  f[2] = c->qm * (E[2] + u[3] * B[1] - u[4] * B[0]);
  PetscCall(VecRestoreArrayRead(X, &x)); PetscCall(VecRestoreArrayRead(c->U, &u)); PetscCall(VecRestoreArray(Vres, &f));
  PetscFunctionReturn(PETSC_SUCCESS);
}
/* Discrete gradient: u' = S(u) grad H(u), S = [[0, I], [-I, qm*Bhat(x)]], Bhat v = v x B (skew) */
static PetscErrorCode DGS(TS ts, PetscReal t, Vec U, Mat S, void *ctx)
{
  Ctx *c = (Ctx *)ctx; const PetscScalar *u; PetscReal E[3], B[3], phi, psi;
  PetscFunctionBeginUser;
  PetscCall(MatSetOption(S, MAT_NEW_NONZERO_ALLOCATION_ERR, PETSC_FALSE));
  PetscCall(VecGetArrayRead(U, &u));
  Fields(c, u, E, B, &phi, &psi);
  PetscCall(VecRestoreArrayRead(U, &u));
  for (PetscInt i = 0; i < 3; i++) {
    PetscCall(MatSetValue(S, i, 3 + i, 1.0, INSERT_VALUES));
    PetscCall(MatSetValue(S, 3 + i, i, -1.0, INSERT_VALUES));
  }
  /* (v x B)_0 = v1 B2 - v2 B1 etc. */
  PetscCall(MatSetValue(S, 3, 4, c->qm * B[2], INSERT_VALUES)); PetscCall(MatSetValue(S, 3, 5, -c->qm * B[1], INSERT_VALUES));
  PetscCall(MatSetValue(S, 4, 5, c->qm * B[0], INSERT_VALUES)); PetscCall(MatSetValue(S, 4, 3, -c->qm * B[2], INSERT_VALUES));
  PetscCall(MatSetValue(S, 5, 3, c->qm * B[1], INSERT_VALUES)); PetscCall(MatSetValue(S, 5, 4, -c->qm * B[0], INSERT_VALUES));
  PetscCall(MatAssemblyBegin(S, MAT_FINAL_ASSEMBLY)); PetscCall(MatAssemblyEnd(S, MAT_FINAL_ASSEMBLY));
  PetscFunctionReturn(PETSC_SUCCESS);
}
static PetscErrorCode DGF(TS ts, PetscReal t, Vec U, PetscScalar *F, void *ctx)
{
  const PetscScalar *u; PetscReal H, P;
  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u)); Invariants((Ctx *)ctx, u, &H, &P); PetscCall(VecRestoreArrayRead(U, &u));
  *F = H;
  PetscFunctionReturn(PETSC_SUCCESS);
}
static PetscErrorCode DGG(TS ts, PetscReal t, Vec U, Vec G, void *ctx)
{
  Ctx *c = (Ctx *)ctx; const PetscScalar *u; PetscScalar *g; PetscReal E[3], B[3], phi, psi;
  PetscFunctionBeginUser;
  PetscCall(VecGetArrayRead(U, &u)); PetscCall(VecGetArray(G, &g));
  Fields(c, u, E, B, &phi, &psi);
  for (PetscInt i = 0; i < 3; i++) { g[i] = -c->qm * E[i]; g[3 + i] = u[3 + i]; }
  PetscCall(VecRestoreArrayRead(U, &u)); PetscCall(VecRestoreArray(G, &g));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* Boris (symmetric drift-kick-drift form): volume preserving, exactly |v|-preserving rotation */
static void BorisStep(Ctx *c, PetscReal *u, PetscReal dt)
{
  PetscReal E[3], B[3], phi, psi, h = 0.5 * dt * c->qm, vm[3], vp[3], t[3], s[3], vprime[3], t2;
  for (int i = 0; i < 3; i++) u[i] += 0.5 * dt * u[3 + i];
  Fields(c, u, E, B, &phi, &psi);
  for (int i = 0; i < 3; i++) { vm[i] = u[3 + i] + h * E[i]; t[i] = h * B[i]; }
  t2 = t[0] * t[0] + t[1] * t[1] + t[2] * t[2];
  for (int i = 0; i < 3; i++) s[i] = 2 * t[i] / (1 + t2);
  vprime[0] = vm[0] + vm[1] * t[2] - vm[2] * t[1];
  vprime[1] = vm[1] + vm[2] * t[0] - vm[0] * t[2];
  vprime[2] = vm[2] + vm[0] * t[1] - vm[1] * t[0];
  vp[0] = vm[0] + vprime[1] * s[2] - vprime[2] * s[1];
  vp[1] = vm[1] + vprime[2] * s[0] - vprime[0] * s[2];
  vp[2] = vm[2] + vprime[0] * s[1] - vprime[1] * s[0];
  for (int i = 0; i < 3; i++) u[3 + i] = vp[i] + h * E[i];
  for (int i = 0; i < 3; i++) u[i] += 0.5 * dt * u[3 + i];
}

static PetscErrorCode Integrate(Ctx *c, const char *method, const PetscReal *u0, PetscReal dt, PetscInt steps, PetscReal *uf)
{
  PetscFunctionBeginUser;
  if (!strcmp(method, "boris")) {
    PetscReal u[6];
    for (int i = 0; i < 6; i++) u[i] = u0[i];
    Record(c, 0, 0, u);
    for (PetscInt n = 1; n <= steps; n++) { BorisStep(c, u, dt); Record(c, n, n * dt, u); }
    for (int i = 0; i < 6; i++) uf[i] = u[i];
  } else {
    TS ts; Vec U; PetscScalar *ua;
    PetscCall(VecCreateSeq(PETSC_COMM_SELF, 6, &U));
    PetscCall(VecGetArray(U, &ua)); for (int i = 0; i < 6; i++) ua[i] = u0[i]; PetscCall(VecRestoreArray(U, &ua));
    c->U = U;
    PetscCall(TSCreate(PETSC_COMM_SELF, &ts));
    PetscCall(TSSetProblemType(ts, TS_NONLINEAR));
    if (!strcmp(method, "rk4")) {
      PetscCall(TSSetType(ts, TSRK)); PetscCall(TSRKSetType(ts, TSRK4));
      PetscCall(TSSetRHSFunction(ts, NULL, RHS, c));
    } else if (!strcmp(method, "bsi")) {
      IS isx, isv;
      PetscCall(ISCreateStride(PETSC_COMM_SELF, 3, 0, 1, &isx)); PetscCall(ISCreateStride(PETSC_COMM_SELF, 3, 3, 1, &isv));
      PetscCall(TSSetType(ts, TSBASICSYMPLECTIC));
      PetscCall(TSRHSSplitSetIS(ts, "position", isx)); PetscCall(TSRHSSplitSetIS(ts, "momentum", isv));
      PetscCall(TSRHSSplitSetRHSFunction(ts, "position", NULL, RHSPos, c));
      PetscCall(TSRHSSplitSetRHSFunction(ts, "momentum", NULL, RHSMom, c));
      PetscCall(TSSetRHSFunction(ts, NULL, RHS, c));
      PetscCall(ISDestroy(&isx)); PetscCall(ISDestroy(&isv));
    } else if (!strcmp(method, "dg")) {
      PetscCall(TSSetType(ts, TSDISCGRAD));
      PetscCall(TSDiscGradSetFormulation(ts, DGS, DGF, DGG, c));
      PetscCall(TSDiscGradSetType(ts, TS_DG_GONZALEZ));
    } else SETERRQ(PETSC_COMM_SELF, PETSC_ERR_ARG_WRONG, "unknown method %s", method);
    PetscCall(TSSetTimeStep(ts, dt)); PetscCall(TSSetMaxSteps(ts, steps)); PetscCall(TSSetMaxTime(ts, 1e30));
    PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
    PetscCall(TSMonitorSet(ts, Monitor, c, NULL));
    PetscCall(TSSetFromOptions(ts));
    PetscCall(TSSolve(ts, U));
    PetscCall(VecGetArray(U, &ua)); for (int i = 0; i < 6; i++) uf[i] = ua[i]; PetscCall(VecRestoreArray(U, &ua));
    PetscCall(TSDestroy(&ts)); PetscCall(VecDestroy(&U));
  }
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscReal Det6(PetscReal A[6][6])
{
  PetscReal d = 1;
  for (int k = 0; k < 6; k++) {
    int p = k;
    for (int i = k + 1; i < 6; i++) if (PetscAbsReal(A[i][k]) > PetscAbsReal(A[p][k])) p = i;
    if (p != k) { for (int j = 0; j < 6; j++) { PetscReal t = A[k][j]; A[k][j] = A[p][j]; A[p][j] = t; } d = -d; }
    if (A[k][k] == 0) return 0;
    d *= A[k][k];
    for (int i = k + 1; i < 6; i++) { PetscReal f = A[i][k] / A[k][k]; for (int j = k; j < 6; j++) A[i][j] -= f * A[k][j]; }
  }
  return d;
}

int main(int argc, char **argv)
{
  Ctx            c;
  char           method[32] = "boris", fname[PETSC_MAX_PATH_LEN] = "out.dat";
  PetscReal      dt = 0.1, u0[6] = {1.0, 0, 0, 0, 1.0, 0.0}, uf[6], eps = 1e-6;
  PetscInt       steps = 100000, n6 = 6, jac = 0;
  PetscLogDouble t0, t1;

  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  c.field = F_UNIFORM; c.qm = 1.0; c.B0 = 1.0; c.kE = 0.1; c.q0 = 1.5; c.kappa = 1.7; c.nout = 100; c.fp = NULL;
  PetscOptionsBegin(PETSC_COMM_WORLD, NULL, "pusher", NULL);
  PetscCall(PetscOptionsEnum("-field", "field", NULL, FieldTypes, (PetscEnum)c.field, (PetscEnum *)&c.field, NULL));
  PetscCall(PetscOptionsString("-method", "boris|rk4|bsi|dg", NULL, method, method, sizeof(method), NULL));
  PetscCall(PetscOptionsString("-o", "output", NULL, fname, fname, sizeof(fname), NULL));
  PetscCall(PetscOptionsReal("-dt", "dt", NULL, dt, &dt, NULL));
  PetscCall(PetscOptionsReal("-qm", "q/m", NULL, c.qm, &c.qm, NULL));
  PetscCall(PetscOptionsReal("-kE", "E well", NULL, c.kE, &c.kE, NULL));
  PetscCall(PetscOptionsInt("-steps", "steps", NULL, steps, &steps, NULL));
  PetscCall(PetscOptionsInt("-nout", "output every", NULL, c.nout, &c.nout, NULL));
  PetscCall(PetscOptionsInt("-jac", "if >0: det of FD Jacobian of the jac-step map (phase-space volume)", NULL, jac, &jac, NULL));
  PetscCall(PetscOptionsReal("-jac_eps", "FD perturbation", NULL, eps, &eps, NULL));
  PetscCall(PetscOptionsRealArray("-u0", "initial state", NULL, u0, &n6, NULL));
  PetscOptionsEnd();
  if (jac > 0) {
    PetscReal J[6][6], up[6], um[6], fp_[6], fm[6];
    for (int j = 0; j < 6; j++) {
      for (int i = 0; i < 6; i++) { up[i] = u0[i]; um[i] = u0[i]; }
      up[j] += eps; um[j] -= eps;
      PetscCall(Integrate(&c, method, up, dt, jac, fp_));
      PetscCall(Integrate(&c, method, um, dt, jac, fm));
      for (int i = 0; i < 6; i++) J[i][j] = (fp_[i] - fm[i]) / (2 * eps);
    }
    PetscCall(PetscPrintf(PETSC_COMM_WORLD, "JACDET method %s field %s nsteps %" PetscInt_FMT " dt %g det-1 %.6e\n", method, FieldTypes[c.field], jac, (double)dt, (double)(Det6(J) - 1.0)));
    PetscCall(PetscFinalize());
    return 0;
  }
  c.fp = fopen(fname, "w");
  PetscCall(PetscTime(&t0));
  PetscCall(Integrate(&c, method, u0, dt, steps, uf));
  PetscCall(PetscTime(&t1));
  fclose(c.fp);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "method %s field %s steps %" PetscInt_FMT " dt %g wall %.3f s\n", method, FieldTypes[c.field], steps, (double)dt, t1 - t0));
  /* machine-readable verification contract line for the agent pipeline */
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "VERIFY {\"code\": \"pusher\", \"method\": \"%s\", \"field\": \"%s\", \"steps\": %" PetscInt_FMT ", \"invariants\": {\"energy_rel_drift_max\": %.3e, \"p_phi_rel_drift_max\": %.3e}, \"finite\": %s}\n", method, FieldTypes[c.field], steps, (double)c.maxdH, (double)c.maxdP, (PetscIsInfOrNanReal(c.maxdH) || c.maxdH > 1e10) ? "false" : "true"));
  PetscCall(PetscFinalize());
  return 0;
}
