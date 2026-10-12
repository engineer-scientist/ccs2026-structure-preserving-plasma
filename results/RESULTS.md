# Phase 2 results: structure-preserving methods + ALCF pilots
_Grok Bot, last updated Mon Oct 12 2026, ~12 AM ET_

All box results come from PETSc 3.26.0 (serial, optimized build) on the Grok Bot box. ALCF runs were all in the
**debug queue, 1–2 nodes**, project DD_plasma_physics, with files only under `~/grokbot-sp/` on each machine.

## 1. Particle pushers (box): `figs/pusher_{es,uniform,solovev}.png`, `figs/pusher_solovev_orbit.png`
Code: `sp/pusher/pusher.c` (C, PETSc TS). One charged particle, 10^6 steps per run.
Fields: (a) ES: B=0, pendulum in x + harmonic in y,z (a separable Hamiltonian); (b) uniform B_z with a radial
electrostatic well, ω_c Δt = 0.2; (c) analytic Solov'ev tokamak equilibrium (R0=B0=1, κ, q0), ω_c Δt ≈ 0.2, passing orbit.
Invariants: H = |v|²/2 + (q/m)φ; canonical toroidal momentum p_φ = R v_φ + (q/m) ψ (exact for (b),(c) by axisymmetry).

| field | method | steps | max \|ΔH/H\| | max \|Δp_φ/p_φ\| | det(J₁₀₀)−1 (phase-volume) |
|---|---|---|---|---|---|
| ES | RK4 (TSRK 4) | 1e6 | 1.0e-2 (secular) | n/a | −3.5e-6 |
| ES | TSBASICSYMPLECTIC (Verlet) | 1e6 | 1.8e-3 (bounded) | n/a | −7.3e-11 |
| ES | Boris | 1e6 | 1.6e-3 (bounded) | n/a | −4.0e-11 |
| ES | TSDISCGRAD (Gonzalez) | 1e6 | **4.9e-13** | n/a | −1.1e-4 |
| uniform B | RK4 | 1e6 | 5.2e-1 | 1.8e-1 | −1.5e-4 |
| uniform B | TSBASICSYMPLECTIC | blows up < 2,000 steps | ∞ | ∞ | +6.3 |
| uniform B | Boris | 1e6 | 4.8e-4 (bounded) | 1.8e-3 (bounded) | −5.8e-11 |
| uniform B | TSDISCGRAD | 1e6 | **6.2e-12** | **1.4e-12** | +6.3e-10 |
| Solov'ev | RK4 | 1e6 | 1.9e-1 (secular) | 1.5e-4 (secular) | −5.2e-5 |
| Solov'ev | TSBASICSYMPLECTIC | blows up < 2,000 steps | ∞ | ∞ | +4.4 |
| Solov'ev | Boris | 1e6 | **1.3e-13** (E=0 ⇒ \|v\| exact) | 4.6e-6 (bounded) | −3.3e-8 |
| Solov'ev | TSDISCGRAD | 1e6 | 4.6e-10 (grows like √t: set by the SNES tolerance) | 2.8e-6 (bounded) | −4.3e-6 |

Takeaways for the talk:
* RK4 drifts secularly in both invariants; at 10^6 steps the energy is off by 19–52 %.
* Boris (volume-preserving) and TSDISCGRAD (energy-preserving) keep both invariants bounded or at the solver tolerance over 10^6 gyro-steps.
* **PETSc's TSBASICSYMPLECTIC assumes a separable H = T(p) + V(q).** With v×B the Hamiltonian is not separable in (x,v), so the method is
  not symplectic there. It blows up within 2,000 steps, and det J − 1 = O(1) shows it. It works for the ES/Vlasov–Poisson case, which is exactly what ex2 uses.
* TSDISCGRAD conserves energy but is not volume-preserving (det − 1 = 1e-4 for the pendulum). This is the energy-versus-symplecticity trade-off.
* Cost (box, 1 core): Boris 0.05 s per 10^6 steps (hand-coded), RK4 0.9 s, TSDISCGRAD 240–355 s (implicit SNES every step).
* Caveat: one particle, one orbit class per field. TSDISCGRAD accuracy is limited by SNES tolerances (rtol 1e-12, atol 1e-14).

## 2. Divergence-free induction (box): `figs/ct_divB_convergence.png`
Code: `sp/ct/ct2d.py` (numpy) and `sp/ct_stag/ct_stag.c` (**PETSc DMSTAG + TS SSP-RK2**, A_z on vertices, B on faces, div B on cells).
2-D ideal induction, periodic box, field-loop advection (uniform flow and rigid rotation), 128², MUSCL-minmod upwind.
* CT: max |∇_h·B| Δx/max|B₀| ≤ 2.2e-16 (numpy) and 3.2e-24 (DMSTAG) for the whole run. A cell-centred B advected with the same scheme reaches **4.9–5.1e-2** within the first steps.
* DMSTAG and numpy CT give the same magnetic energy (1.402977e-7 → 1.273035e-7) to 7 digits, an independent cross-check.
* Smooth-field convergence (N = 32…256): observed order **p = 2.00, 2.00, 2.00**. With a smooth field the CT and baseline errors are the same, so CT pays no accuracy penalty.
* Caveat: 2-D kinematic induction only (prescribed velocity, no Lorentz-force feedback). The DMSTAG version has been run serially only.

## 3. ALCF pilots: PETSc ex2 (DMSwarm Vlasov–Poisson, Landau damping, TSBASICSYMPLECTIC)
Setup: 1D-1V, k = 0.5, α = 0.01, x ∈ [0, 4π] periodic, v ∈ [−6, 6], P1 FE Poisson with GAMG. Weak scaling at fixed work per rank:
**256 x-cells × 200 v-cells = 51,200 particles per MPI rank** (the x-mesh is refined as ranks are added). Option files: `alcf/scripts/`.
All debug queue, 1–2 nodes. Raw per-job numbers: `alcf/*_summary.json`.

### Did the ALCF runs use GPUs?
**Mostly no.** Be explicit about this on the slides.
* All physics and weak-scaling runs used the **CPU code path**: MPI ranks on CPU cores (Polaris: 32 ranks/node on the EPYC; Aurora: 96 ranks/node on the Xeon Max).
  The Aurora PETSc 3.26 build I made is CPU-only. On Polaris the PETSc build has CUDA, but ex2 uses it only when CUDA vector/matrix types are requested.
* GPU trial on Polaris (job 7738042, 4 ranks, 1 per A100, `-dm_vec_type cuda -dm_mat_type aijcusparse -use_gpu_aware_mpi 0`): runs correctly,
  with total energy identical to the CPU run (ΔW/W = 2.326e-6). The Poisson KSP/SNES ran 100 % on the GPU, but the run is **1.35× slower** (TSStep 34.6 s vs 25.5 s).
  In ex2 only the field solve is offloaded, and it is about 1 % of the work. Particle push, location and migration stay on the CPU.
* GPU trial on Aurora (module petsc/3.25.2-sycl, Kokkos/SYCL, 1 rank because that module fails at >1 rank): the first attempt (job 8921944) failed until `ZES_ENABLE_SYSMAN=1` was set. With it (job 8922008, 1 rank on one PVC tile,
  `-dm_vec_type kokkos -dm_mat_type aijkokkos`), the 51,200-particle run completes, and its final E-field and energy match the CPU run of the same module
  to 4e-14 relative. But TSStep takes **169.6 s versus 14.5 s on one CPU core (11.7× slower)**: the tiny 256-cell Poisson solve is dominated by launch and transfer
  overhead, and particles stay on the host. A larger-mesh GPU run (4,096 cells) was cancelled to save node-hours because the conclusion was already clear.
* Bottom line: getting a real GPU speedup needs the DMSwarm particle kernels (push, deposit, locate) on the device. ex2 as written does not do that.

### Physics/conservation (600 steps, t = 30, Δt = 0.05): `figs/ex2_landau_energy.png`
| machine | ranks | particles | fitted γ (theory −0.1533) | max \|ΔW/W₀\| | max \|ΔW\| / peak field energy |
|---|---|---|---|---|---|
| Polaris | 1 | 51,200 | −0.1613 | 1.37e-5 | 6.8 % |
| Polaris | 4 / 8 / 16 | 0.2 M / 0.4 M / 0.8 M | −0.1612 (all) | 1.30e-5 | 6.5 % |
| Aurora | 1 | 51,200 | −0.1613 | 1.37e-5 | 6.8 % |
| Aurora | 12 | 614,400 | −0.1612 | 1.29e-5 | 6.5 % |
* Aurora np=1 (oneAPI) and Polaris np=1 (GNU) match to every printed digit. Results do not depend on rank count.
* γ is about 5 % steeper than linear theory (envelope fit, t ∈ [0.5, 15.5], 256-cell mesh). On the box with 4,000 particles, −0.156.
* Total energy is **not** exactly conserved: the drift is about 1e-5 of the total but 6–7 % of the peak field energy. Verlet with an FE-projected field is
  symplectic in the particle variables only for a frozen field. Energy-conserving PIC (implicit/discrete-gradient coupling) is the natural next step.
  With Δt = 0.01 (scaling runs) the drift falls to 2.2e-6, about 5× smaller for 5× smaller Δt, roughly first order. Treat that as a hint, not a measured order.

### Weak-scaling slice (100 steps, Δt = 0.01): `figs/ex2_weak_scaling_pilot.png`
| machine | ranks | nodes | TSStep s/step | efficiency T(1)/T(P) | Poisson (SNES) s/step | DMPlexLocatePoints max/min across ranks |
|---|---|---|---|---|---|---|
| Polaris | 1 | 1 | 0.170 | 1.00 | 0.0023 | 1.0 |
| Polaris | 4 / 8 / 16 / 32 | 1 | 0.255 / 0.350 / 0.536 / 0.941 | 0.66 / 0.49 / 0.32 / 0.18 | 0.005–0.013 | – |
| Polaris | 32 (16 per node) | 2 | 0.892 | 0.19 | 0.021 | – |
| Polaris | 64 | 2 | 1.687 | 0.10 | 0.027 | – |
| Aurora | 1 | 1 | 0.128 | 1.00 | 0.0021 | 1.0 |
| Aurora | 12 / 24 / 48 / 96 | 1 | 0.339 / 0.571 / 1.306 / 2.305 | 0.38 / 0.22 / 0.098 / 0.056 | 0.009–0.021 | 3.1 (12), 18.6 (96) |
| Aurora | 192 | 2 | 4.371 | 0.029 | 0.040 | 36.1 |
| Aurora | 96 (48 per node) | 2 | 2.352 | 0.054 | 0.024 | 18.8 |

Findings (each a good talk point):
1. **The field solver scales fine. Particle location does not.** SNES/KSP is 1–2 % of a step at every size. `DMSwarmMigrate` → `DMPlexLocatePoints`
   accounts for nearly all of the time. Its max/min across ranks rises from 1 to 18.6 (96 ranks) and 36 (192 ranks), so a few ranks do most of the
   point location while the others wait. MPI message volume in the migration is small. Hypothesis: the rank(s) holding the periodic-wrap cell get a
   bounding box / hash grid that spans the whole domain and fall back to slow search. That needs one small profiling job to confirm.
2. Both machines follow the **same efficiency-vs-P curve**, so the loss is algorithmic, not a hardware or network effect. Spreading 32 ranks over
   2 nodes buys only 5 %.
3. **Silent particle loss (a correctness bug worth reporting):** with Δt = 0.05, runs with ≥ 48 ranks lost 1 %–35 % of particles in the first step
   (Polaris 64 ranks: 2.98 M of 3.28 M left), and the solver then failed. The cause is that the |v|max Δt = 0.3 drift exceeds the subdomain width 4π/P,
   and DMSwarm migration only reaches neighbouring ranks. With Δt = 0.01 the particle count is exactly conserved up to 192 ranks. This is a natural
   check for the VERIFY{} contract: total particle count and charge as invariants.
4. In 1D with a fixed box, weak scaling is capped near P_max ≈ L/(|v|max Δt). A larger study should grow the box with P (L = 4π·P/P₀, k = 0.5 unchanged).

## 4. VERIFY{} conservation contract (agent pipeline prototype)
Repo copy `/workspace/research/petsc_mcp_servers_verify`, local branch `grokbot/verify-contract`, **not pushed** (patches in `code_bundle/`).
* `verify_contract.py`: generated codes must print `VERIFY {"code":…, "invariants":{"energy_rel_drift_max":…, "p_phi_rel_drift_max":…, "divB_max":…}, "finite":…}`. `check()` compares the result against per-invariant tolerances.
* `petsc_compile_run_mcp_server.py`: new `run_and_verify` tool. The code-generator prompt now requires the VERIFY line and loops until verification passes.
* Offline test (`test_verify_contract.py`, real pusher runs, Solov'ev): Boris **PASS**, TSDISCGRAD **PASS**, RK4 **FAIL** (ΔH = 1.5e-2 > 1e-6),
  TSBASICSYMPLECTIC in a B field **FAIL** (non-finite). `ct_stag` also emits a VERIFY line (div B = 3e-24, PASS).
* Caveat: the end-to-end LLM loop was **not** run (no LLM backend is configured outside ANL). Only the deterministic contract and tool code were tested.

## ALCF usage (walltime × nodes, from `qstat -x`)
* Polaris: 7737974 (1 node, 30:14), 7738002 (2 nodes, 6:53), 7738042 (2 nodes, 16:10) → **1.27 node-hours**
* Aurora: 8921827 (1, 3:41, cancelled), 8921836 (1, 3:33), 8921861 (2, 0:23), 8921869 (1, 1:02), 8921873 (1, 0:16), 8921878 (1, 11:32),
  8921888 (2, 0:25), 8921903 (1, 3:34), 8921916 (1, 2:49), 8921920 (1, 14:54), 8921924 (2, 22:23), 8921944 (1, 7:51; 8922008 (1, 3:59, GPU)) → **1.66 node-hours**
* Debug queue only, never more than 2 nodes. Builds were done on login nodes in HOME (~/grokbot-sp).

## Proposal for a larger weak-scaling run (needs your approval; nothing submitted)
Recommendation: **do not spend a large allocation until the point-location imbalance is fixed**. Right now a big run would only re-measure the same pathology.
1. Step A (within current approval, ≤ 0.5 node-h): one 1-node Aurora debug job with per-rank `-log_view` detail and two variants
   (`-dm_plex_hash_location 0` on x; a non-periodic-localized coordinate test) to confirm the hot rank. Then patch ex2 or the options locally (no push).
2. Step B (**needs approval**): grow-the-box weak scaling (L = 4π·P/48, 256 cells and 51,200 particles per rank, Δt = 0.05, 100 steps, about 5 min per run):
   * Aurora `debug-scaling` queue: 1, 2, 4, 8, 16 nodes at 96 ranks/node (96 → 1,536 ranks), one job each, ≤ 20 min walltime:
     about 31 nodes × 0.33 h ≈ **10 node-hours**.
   * Polaris `debug-scaling` queue: 1, 2, 4, 8 nodes at 32 ranks/node (+ the same at 4 ranks/node with CUDA vectors for comparison):
     15 nodes × 0.33 h × 2 ≈ **10 node-hours**.
   * Total ≈ **20 node-hours** of the DD_plasma_physics allocation. Please check the current per-queue node and walltime limits before approving.
   * Expected value: (a) one weak-scaling curve per machine out to 8–16 nodes for the talk or paper, with the box grown so the physics stays valid;
     (b) a before/after demonstration of the locate-imbalance fix, if Step A finds one; (c) particle count, charge and energy invariants checked through the
     VERIFY{} contract at every scale. Without Step A the expected value is low, since the efficiency would just keep falling as in the current curve.
3. GPU-resident particles (DMSwarm on device) are a development task, not a run, so they are not budgeted here.

## Honest caveats / issues found
* Aurora module `petsc/3.25.2-sycl`: ex2 runs at 1 rank, but every multi-rank run aborted with `MPI_Abort(…, 63)` and no PETSc message (not root-caused). Multi-rank Aurora results use my CPU-only PETSc 3.26.0 build. On Aurora the GPU was exercised only at 1 rank via the module (Kokkos/SYCL, needs `ZES_ENABLE_SYSMAN=1`).
* My first Poisson settings (GMRES, ksp atol 1e-14) and the v2 settings (CG+GAMG, atol 1e-12) both failed in runs that had lost particles. With Δt = 0.01, GMRES+GAMG at rtol 1e-8 / atol 1e-9 ran cleanly up to 192 ranks.
* The weak-scaling runs are only 100 steps (t = 1) and cover only 1D-1V. The 600-step physics runs go up to 16 ranks (Polaris) and 12 ranks (Aurora).
* Earlier failed or cancelled jobs are listed under ALCF usage. Every result reported here comes from a completed run.
* An ex2 GPU run means only PETSc CUDA vectors/matrices for the mesh DM. DMSwarm particle kernels in ex2 run on the CPU, so GPU speedups are not expected from this code as-is.

## Phase 3 (Oct 12, 2026): locate imbalance diagnosed and fixed, then weak scaling to 16 nodes
Approved 12:40 AM ET Oct 12. Running log: `phase3/NOTES.md`. Scripts: `phase3/scripts/` (copies in `~/grokbot-sp/p3/` on both machines).
Figure: `figs/weak_scaling_v2.png`. Data: `weak_scaling_v2.csv` (combined), `aurora_p3diag_8922161_reduced.csv`, `aurora_p3ws_8922217_reduced.csv`, `polaris_p3ws_7738212_reduced.csv`.

### Step A: why one rank was 18–41× slower in DMPlexLocatePoints
* **Cause (confirmed three ways).** `DMPlexCreateGridHash()` (PETSc 3.26 `plexgeometry.c`) sizes each rank's point-location hash box from *all* of that rank's local coordinates.
  On a periodic 1D mesh, the rank that owns the wrap-around cell (the last rank) also holds the vertex at x = 0, so its hash box spans the
  **whole domain [0, L]** while it gets only ~0.8 × 256 ≈ 204 bins. Its 256 cells then fall into one or two bins, and every particle that left its
  previous cell costs a linear scan over ~256 cells instead of ~1.
  1. `-dm_plex_print_locate` on the box (4 ranks): ranks 0–2 report hash boxes equal to their own subdomain; rank 3 reports `(0 -- 12.566)`.
  2. Per-rank `-log_view :…:ascii_info_detail` on Aurora (job 8922161, 96 ranks): DMPlexLocatePoints time **131 s on rank 95 vs a median of 7.4 s** (18.6×);
     flops on rank 95 are 19× those on any other rank. Every other rank waits for it in the next collective, so the whole step runs at rank 95's speed.
  3. The imbalance follows the formula: with 256·P/204 cells per bin on the wrap rank, it grows with P until a single bin holds all 256 cells,
     which matches the saturation at **≈ 41×** seen at 384–1,536 ranks.
  It is not the initial particle distribution (particles are created per local cell and stay balanced), not periodic coordinate localization itself,
  and not migration volume (message sizes are small and balanced). Brute-force location (`-dm_plex_hash_location 0`) makes it worse: 5.3 s per step at 96 ranks.
* **Fix.** A small local patch to ex2 (`phase3/patch/ex2_gb_wrap_hash.patch`, 27 lines, not pushed): new option `-gb_wrap_hash_factor f`. After the mesh
  is created, any rank whose local bounding box is more than 4× wider than its owned cells (only the wrap rank) sets `-dm_plex_hash_box_faces = f·width/h`
  in its own options database, so its bins are about one cell wide. All other ranks keep the default hash. We used f = 2.
  The upstream fix would be for `DMPlexCreateGridHash()` to size the box from the localized coordinates of owned cells only. That is worth reporting to PETSc.
* **Rejected first fix.** The global option `-dm_plex_hash_box_faces 2·256·P` works at 96 ranks (8.2× faster), but the hash build costs grow faster
  than linearly with the number of bins (box test: 4k bins 0.3 s, 65k 8 s, 262k 145 s). At 1,536 ranks (786k bins on *every* rank) setup never finished,
  and job 8922177 was cancelled after 10 minutes (2.7 node-hours lost). The patch gives fine bins only to the one rank that needs them.
* Effect at 96 ranks (one Aurora node): per-step time **2.29 s → 0.28–0.30 s (7.7–8.2× faster)**. Locate work imbalance (max/min flops) **19× → 1.07–1.17×**.

### Step B: grow-the-box weak scaling with the fix
Setup: 1D-1V Landau damping, k = 0.5, 256 x-cells × 200 v-cells = **51,200 particles per rank**, Δt = 0.05, 100 steps ("before": 20 steps).
**Domain L = π·P** (4π·P/4, a width of π per rank). This differs from the proposal's L = 4π·P/48: that rule gives a subdomain width of 0.26,
below |v|max·Δt = 0.3, so particles would skip a rank and be lost. With width π, every P that is a multiple of 4 holds a whole number of k = 0.5 wavelengths.
Efficiency below is T(1 rank)/T(P) from the TSStep time per step (slowest rank). The 1-node-relative value is in brackets.

| Aurora (96 ranks/node) | 1 rank | 96 (1 node) | 192 (2) | 384 (4) | 768 (8) | 1,536 (16) |
|---|---|---|---|---|---|---|
| s/step after fix | 0.147 | 0.299 | 0.330 | 0.342 | 0.363 | 0.400 |
| efficiency after fix [vs 1 node] | 1.00 | 0.49 [1.00] | 0.45 [0.91] | 0.43 [0.87] | 0.41 [0.82] | 0.37 [0.75] |
| s/step before fix | 0.151 | 2.29 | 4.32 | 4.93 | 4.85 | 4.88 |
| efficiency before fix | 1.00 | 0.066 | 0.035 | 0.031 | 0.031 | 0.031 |
| Locate max/min across ranks, before → after | – | 18.2 → 1.2 | 35.6 → 1.2 | 41.5 → 1.2 | 41.5 → 1.4 | 41.4 → 1.3 |
| particles (all conserved exactly) | 51,200 | 4.9 M | 9.8 M | 19.7 M | 39.3 M | **78.6 M** |

| Polaris (job 7738212, 8 nodes, debug-scaling) | 1 rank | 32 (1 node) | 64 (2) | 128 (4) | 256 (8) |
|---|---|---|---|---|---|
| s/step after fix, 32 CPU ranks/node | 0.192 | 0.269 | 0.291 | 0.314 | 0.337 |
| efficiency after fix [vs 1 node] | 1.00 | 0.71 [1.00] | 0.66 [0.92] | 0.61 [0.86] | 0.57 [0.80] |
| s/step before fix | – | 0.949 | 1.721 | 2.798 | 4.958 |
| efficiency before fix [vs 1 node] | – | 0.20 [1.00] | 0.11 [0.55] | 0.07 [0.34] | 0.04 [0.19] |
| Locate max/min before → after | – | 6.9 → 1.2 | 13.2 → 1.1 | 22.7 → 1.2 | 41.4 → 1.1 |
| particles (all conserved exactly) | 51,200 | 1.6 M | 3.3 M | 6.6 M | **13.1 M** |

| Polaris, 4 ranks/node (1 per A100) | 4 (1 node) | 8 (2) | 16 (4) | 32 (8) |
|---|---|---|---|---|
| CPU: s/step [vs 1 node] | 0.278 [1.00] | 0.230 [1.21] | 0.234 [1.19] | 0.248 [1.12] |
| CUDA vec/mat: s/step [vs 1 node] | 0.357 [1.00] | 0.336 [1.06] | 0.340 [1.05] | 0.401 [0.89] |
| GPU slowdown vs CPU at the same ranks | 1.29× | 1.46× | 1.46× | 1.62× |

* Job 7738212 includes all of the earlier partial job 7738293 (1–2 nodes) and agrees with it to within 5%. It replaces 7738293 in the table, the figure and the CSV.
  (At 4 ranks/node, efficiency vs 1 node is above 1: the 1-node run takes 0.278 s/step, one of four ranks has slower locate times, and that does not repeat at larger counts.)
* GPU (CUDA vectors/matrices, 1 rank per A100) is still **1.3–1.6× slower** than the same number of CPU ranks (the gap grows with node count), because only the field solve is offloaded (as in Phase 2).

* **Speed-up from the fix at 16 Aurora nodes: 12.2× per step.** Across nodes the fixed code holds 75% efficiency from 1 to 16 nodes on Aurora and 80% from 1 to 8 nodes on Polaris (12.2× and 14.7× faster than before the fix at 16 Aurora and 8 Polaris nodes).
  Before the fix, efficiency was 3% and was set entirely by the one wrap rank.
* **What is left.** (i) A drop within one node (1 rank → full node: 0.49 on Aurora). The point-location time per rank roughly doubles (14 s → 25 s per 100 steps)
  even though the work is balanced, which points to shared memory bandwidth and cache, not imbalance. (ii) A slow rise from 1 to 16 nodes, mostly the Poisson solve
  (SNESSolve 0.024 → 0.073 s/step, GAMG on a 393k-cell 1D mesh) and mass-matrix assembly. Location is now flat with P.
* **Invariants (the VERIFY{} contract) at every scale.** Global particle count and total weight (charge) are identical at the first and last step in every run,
  1 to 1,536 ranks on Aurora and 1 to 256 ranks on Polaris (CPU and GPU), with no particle loss. The relative total-energy drift is the same at every P: **1.14e-5** after 100 steps (Δt = 0.05), and 1.12e-5 for the 20-step runs.
  So the parallel decomposition does not change the physics. (The 1-rank run, L = π, has a different box and gives 1.1e-7.)

### ALCF usage, Phase 3 (walltime × nodes from `qstat -xf`)
* Aurora: 8922147 (1 node, 0:03, failed: bad symlink) ≈ 0.00; 8922161 diagnostic (1 node, 8:27) = 0.14; **8922177 (16 nodes, 10:12, stalled in hash setup, cancelled) = 2.72 (wasted)**;
  8922217 sweep (16 nodes, 19:33) = 5.21 → **Aurora 8.08 node-hours**
* Polaris: 7738204 (cancelled before it started) = 0; 7738288 (2 nodes, 0:02, failed: wrong MPI library at link time) ≈ 0.00; 7738293 (2 nodes, 6:54) = 0.23; 7738212 (8 nodes, 17:09) = 2.29 → **Polaris 2.52**
* **Phase 3 total (final): 10.60 node-hours** of the 25 cap, of which 2.72 were wasted on 8922177. No jobs remain.

### Caveats
* The fix is a local ex2 patch, not a PETSc change. Other DMSwarm codes on periodic Plex meshes in 2D and 3D likely have the same wrap-rank hash problem,
  and the patch as written only handles 1D (`cdim == 1`).
* "Before" points in panel (a) come from two setups: Phase 2 (fixed box L = 4π, Δt = 0.01, mesh refined with P) and Phase 3 (grow box, Δt = 0.05, 20 steps).
  Both show the same pathology.
* Runs are short (100 steps, t = 5). Times exclude setup. At 1,536 ranks the wrap rank spends roughly 20 s once, at startup, building its fine hash.
* The ex2 "GPU" runs offload only the field solve. Particles stay on the CPU, as in Phase 2.
