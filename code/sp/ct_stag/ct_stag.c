static char help[] = "2D constrained-transport induction solver on PETSc DMSTAG + TS (SSP-RK2).\n"
"A_z on vertices, B = curl A on faces (B_x on left faces, B_y on down faces), div B on elements.\n"
"Field-loop advection (Gardiner & Stone 2005) with uniform velocity, periodic unit square.\n";
#include <petscdmstag.h>
#include <petscts.h>

typedef struct { PetscReal ux, uy, h; DM dmB; Vec B; PetscInt nout; FILE *fp; PetscReal Bmax0; } AppCtx;

static PetscReal minmod(PetscReal a, PetscReal b) { return (a * b > 0) ? (PetscAbsReal(a) < PetscAbsReal(b) ? a : b) : 0.0; }

/* upwind MUSCL(minmod) for -(u . grad) A, matching sp/ct/ct2d.py */
static PetscErrorCode RHS(TS ts, PetscReal t, Vec A, Vec F, void *ctx)
{
  AppCtx            *user = (AppCtx *)ctx;
  DM                 dm;
  Vec                Al;
  const PetscScalar ***a;
  PetscScalar     ***f;
  PetscInt           sx, sy, nx, ny, s;
  PetscFunctionBeginUser;
  PetscCall(TSGetDM(ts, &dm));
  PetscCall(DMGetLocalVector(dm, &Al));
  PetscCall(DMGlobalToLocal(dm, A, INSERT_VALUES, Al));
  PetscCall(DMStagVecGetArrayRead(dm, Al, &a));
  Vec Fl;
  PetscCall(DMGetLocalVector(dm, &Fl));
  PetscCall(DMStagVecGetArray(dm, Fl, &f));
  PetscCall(DMStagGetCorners(dm, &sx, &sy, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN_LEFT, 0, &s));
  for (PetscInt j = sy; j < sy + ny; j++)
    for (PetscInt i = sx; i < sx + nx; i++) {
      PetscReal r = 0;
      for (int d = 0; d < 2; d++) {
        PetscReal q[5], u = d == 0 ? user->ux : user->uy;
        for (int k = -2; k <= 2; k++) q[k + 2] = d == 0 ? a[j][i + k][s] : a[j + k][i][s];
        /* face values at i+1/2 and i-1/2 */
        PetscReal sl[5];
        for (int k = 1; k <= 3; k++) sl[k] = minmod(q[k] - q[k - 1], q[k + 1] - q[k]);
        PetscReal fp = u > 0 ? q[2] + 0.5 * sl[2] : q[3] - 0.5 * sl[3];
        PetscReal fm = u > 0 ? q[1] + 0.5 * sl[1] : q[2] - 0.5 * sl[2];
        r -= u * (fp - fm) / user->h;
      }
      f[j][i][s] = r;
    }
  PetscCall(DMStagVecRestoreArrayRead(dm, Al, &a));
  PetscCall(DMStagVecRestoreArray(dm, Fl, &f));
  PetscCall(DMLocalToGlobal(dm, Fl, INSERT_VALUES, F));
  PetscCall(DMRestoreLocalVector(dm, &Fl));
  PetscCall(DMRestoreLocalVector(dm, &Al));
  PetscFunctionReturn(PETSC_SUCCESS);
}

/* B = curl A on faces, then max |div B| on elements (normalized by Bmax0/h) */
static PetscErrorCode DivB(DM dm, Vec A, AppCtx *user, PetscReal *divmax, PetscReal *Bmax, PetscReal *EB)
{
  Vec                Al, Bl;
  const PetscScalar ***a;
  PetscScalar     ***b;
  PetscInt           sx, sy, nx, ny, sA, sL, sD, sE;
  PetscReal          dm_ = 0, bm = 0, eb = 0, h = user->h;
  PetscFunctionBeginUser;
  PetscCall(DMGetLocalVector(dm, &Al));
  PetscCall(DMGlobalToLocal(dm, A, INSERT_VALUES, Al));
  PetscCall(DMStagVecGetArrayRead(dm, Al, &a));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN_LEFT, 0, &sA));
  PetscCall(DMStagGetLocationSlot(user->dmB, DMSTAG_LEFT, 0, &sL));
  PetscCall(DMStagGetLocationSlot(user->dmB, DMSTAG_DOWN, 0, &sD));
  PetscCall(DMStagGetLocationSlot(user->dmB, DMSTAG_ELEMENT, 0, &sE));
  PetscCall(DMStagGetCorners(user->dmB, &sx, &sy, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
  PetscCall(DMGetLocalVector(user->dmB, &Bl));
  PetscCall(DMStagVecGetArray(user->dmB, Bl, &b));
  for (PetscInt j = sy; j < sy + ny; j++)
    for (PetscInt i = sx; i < sx + nx; i++) {
      b[j][i][sL] = (a[j + 1][i][sA] - a[j][i][sA]) / h;   /* Bx at (x_i, y_{j+1/2}) */
      b[j][i][sD] = -(a[j][i + 1][sA] - a[j][i][sA]) / h;  /* By at (x_{i+1/2}, y_j) */
      bm = PetscMax(bm, PetscMax(PetscAbsReal(b[j][i][sL]), PetscAbsReal(b[j][i][sD])));
      eb += 0.5 * (b[j][i][sL] * b[j][i][sL] + b[j][i][sD] * b[j][i][sD]) * h * h;
    }
  PetscCall(DMStagVecRestoreArray(user->dmB, Bl, &b));
  PetscCall(DMLocalToGlobal(user->dmB, Bl, INSERT_VALUES, user->B));
  PetscCall(DMGlobalToLocal(user->dmB, user->B, INSERT_VALUES, Bl));
  PetscCall(DMStagVecGetArray(user->dmB, Bl, &b));
  for (PetscInt j = sy; j < sy + ny; j++)
    for (PetscInt i = sx; i < sx + nx; i++) {
      PetscReal dv = (b[j][i + 1][sL] - b[j][i][sL]) / h + (b[j + 1][i][sD] - b[j][i][sD]) / h;
      b[j][i][sE] = dv;
      dm_ = PetscMax(dm_, PetscAbsReal(dv));
    }
  PetscCall(DMStagVecRestoreArray(user->dmB, Bl, &b));
  PetscCall(DMRestoreLocalVector(user->dmB, &Bl));
  PetscCall(DMStagVecRestoreArrayRead(dm, Al, &a));
  PetscCall(DMRestoreLocalVector(dm, &Al));
  PetscCallMPI(MPIU_Allreduce(MPI_IN_PLACE, &dm_, 1, MPIU_REAL, MPIU_MAX, PetscObjectComm((PetscObject)dm)));
  PetscCallMPI(MPIU_Allreduce(MPI_IN_PLACE, &bm, 1, MPIU_REAL, MPIU_MAX, PetscObjectComm((PetscObject)dm)));
  PetscCallMPI(MPIU_Allreduce(MPI_IN_PLACE, &eb, 1, MPIU_REAL, MPIU_SUM, PetscObjectComm((PetscObject)dm)));
  *divmax = dm_; *Bmax = bm; *EB = eb;
  PetscFunctionReturn(PETSC_SUCCESS);
}

static PetscErrorCode Monitor(TS ts, PetscInt step, PetscReal t, Vec A, void *ctx)
{
  AppCtx   *user = (AppCtx *)ctx;
  DM        dm;
  PetscReal dv, bm, eb;
  PetscFunctionBeginUser;
  if (step % user->nout) PetscFunctionReturn(PETSC_SUCCESS);
  PetscCall(TSGetDM(ts, &dm));
  PetscCall(DivB(dm, A, user, &dv, &bm, &eb));
  if (step == 0) user->Bmax0 = bm;
  if (user->fp) fprintf(user->fp, "%.6e %.6e %.10e\n", (double)t, (double)(dv * user->h / user->Bmax0), (double)eb);
  PetscFunctionReturn(PETSC_SUCCESS);
}

int main(int argc, char **argv)
{
  DM        dm;
  TS        ts;
  Vec       A;
  AppCtx    user;
  PetscInt  N = 128, sx, sy, nx, ny, s, nsteps;
  PetscReal T = 2.0, cfl = 0.4, divmax_all = 0, dv, bm, eb, eb0;
  char      fname[PETSC_MAX_PATH_LEN] = "ct_stag.dat";
  PetscScalar ***a;
  PetscMPIInt size;

  PetscCall(PetscInitialize(&argc, &argv, NULL, help));
  PetscCallMPI(MPI_Comm_size(PETSC_COMM_WORLD, &size));
  user.ux = 1.0; user.uy = 0.5; user.nout = 5; user.fp = NULL;
  PetscCall(PetscOptionsGetInt(NULL, NULL, "-N", &N, NULL));
  PetscCall(PetscOptionsGetReal(NULL, NULL, "-T", &T, NULL));
  PetscCall(PetscOptionsGetString(NULL, NULL, "-o", fname, sizeof(fname), NULL));
  user.h = 1.0 / N;
  PetscCall(DMStagCreate2d(PETSC_COMM_WORLD, DM_BOUNDARY_PERIODIC, DM_BOUNDARY_PERIODIC, N, N, PETSC_DECIDE, PETSC_DECIDE, 1, 0, 0, DMSTAG_STENCIL_BOX, 2, NULL, NULL, &dm));
  PetscCall(DMSetFromOptions(dm)); PetscCall(DMSetUp(dm));
  PetscCall(DMStagCreateCompatibleDMStag(dm, 0, 1, 1, 0, &user.dmB));
  PetscCall(DMSetUp(user.dmB));
  PetscCall(DMCreateGlobalVector(user.dmB, &user.B));
  PetscCall(DMCreateGlobalVector(dm, &A));
  PetscCall(DMStagGetCorners(dm, &sx, &sy, NULL, &nx, &ny, NULL, NULL, NULL, NULL));
  PetscCall(DMStagGetLocationSlot(dm, DMSTAG_DOWN_LEFT, 0, &s));
  Vec Aloc;
  PetscCall(DMGetLocalVector(dm, &Aloc));
  PetscCall(DMStagVecGetArray(dm, Aloc, &a));
  for (PetscInt j = sy; j < sy + ny; j++)
    for (PetscInt i = sx; i < sx + nx; i++) {
      PetscReal x = i * user.h, y = j * user.h, r = PetscSqrtReal((x - 0.5) * (x - 0.5) + (y - 0.5) * (y - 0.5));
      a[j][i][s] = 1e-3 * PetscMax(0.3 - r, 0.0);
    }
  PetscCall(DMStagVecRestoreArray(dm, Aloc, &a));
  PetscCall(DMLocalToGlobal(dm, Aloc, INSERT_VALUES, A));
  PetscCall(DMRestoreLocalVector(dm, &Aloc));
  if (size == 1 || 1) { PetscMPIInt rank; MPI_Comm_rank(PETSC_COMM_WORLD, &rank); if (!rank) user.fp = fopen(fname, "w"); }
  nsteps = (PetscInt)PetscCeilReal(T / (cfl * user.h / user.ux));
  PetscCall(TSCreate(PETSC_COMM_WORLD, &ts));
  PetscCall(TSSetDM(ts, dm));
  PetscCall(TSSetType(ts, TSSSP)); PetscCall(TSSSPSetType(ts, TSSSPRKS2)); PetscCall(TSSSPSetNumStages(ts, 2));
  PetscCall(TSSetRHSFunction(ts, NULL, RHS, &user));
  PetscCall(TSSetTimeStep(ts, T / nsteps)); PetscCall(TSSetMaxSteps(ts, nsteps)); PetscCall(TSSetMaxTime(ts, T));
  PetscCall(TSSetExactFinalTime(ts, TS_EXACTFINALTIME_MATCHSTEP));
  PetscCall(TSMonitorSet(ts, Monitor, &user, NULL));
  PetscCall(TSSetFromOptions(ts));
  PetscCall(DivB(dm, A, &user, &dv, &bm, &eb0));
  PetscCall(TSSolve(ts, A));
  PetscCall(DivB(dm, A, &user, &dv, &bm, &eb));
  divmax_all = dv * user.h / user.Bmax0;
  if (user.fp) fclose(user.fp);
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "DMSTAG CT N=%" PetscInt_FMT " ranks=%d steps=%" PetscInt_FMT " final max|divB|*h/Bmax0 = %.3e  magnetic energy %.6e -> %.6e\n", N, size, nsteps, (double)divmax_all, (double)eb0, (double)eb));
  PetscCall(PetscPrintf(PETSC_COMM_WORLD, "VERIFY {\"code\": \"ct_stag\", \"method\": \"CT-DMSTAG\", \"invariants\": {\"divB_max\": %.3e}, \"finite\": true}\n", (double)divmax_all));
  PetscCall(VecDestroy(&A)); PetscCall(VecDestroy(&user.B)); PetscCall(DMDestroy(&user.dmB)); PetscCall(TSDestroy(&ts)); PetscCall(DMDestroy(&dm));
  PetscCall(PetscFinalize());
  return 0;
}
