# Research Plan — Structure-Preserving Plasma Simulation in PETSc + Multi-Agent Generation + ALCF Weak Scaling

Owner: Sarthak Sharma (UB; w/ M. Knepley). Prepared 2026-10-11 (ET) by Grok Bot, Phase 1 (read-only survey).
Target: results before ALCF DD allocation `DD_plasma_physics` ends **2026-11-07**; then a 12-min virtual talk at CCS 2026.

---------------------------------------------------------------------------------------------------
## 0. What exists today (grounded in what was actually read)

### 0.1 `petsc_mcp_servers` (gitlab.com/petsc/petsc_mcp_servers, cloned at /workspace/research/petsc_mcp_servers, HEAD 21ccac2)
FastMCP servers, each an LLM-backed (Claude Agent SDK) or deterministic tool server:
| Server | Port | What it really does |
|---|---|---|
| `orchestrator_mcp_server` | 8086 | Inner Claude told to call, in fixed order: PDE modeling -> NA -> code generator -> compile-run. `cntlimit=35` message cap. No verification step of its own. |
| `pde_modeling_mcp_server` | 8084 | LLM returns PDE in HTML/LaTeX/UFL + "time-dependent" flag + name. |
| `na_mcp_server` | 8085 | LLM returns JSON {grid, discretization, integrator}; deterministic maps grid/disc -> DMDA/DMSTAG/DMPLEX/DMSWARM and problem -> TS/SNES/KSP. **No notion of structure preservation, conservation laws, or symplecticity.** |
| `petsc_claude_code_generator_mcp_server` | 8083 | LLM writes C, uses compile-run + docs + RAG; then server re-checks it *compiles* (not correctness). |
| `petsc_compile_run_mcp_server` | — | bash-like sandbox: create file, make, run_executable (20 s default timeout, local mpiexec). No PBS/HPC submission. |
| `petsc_documentation_mcp_server`, `petsc_rag_mcp_server` | — | PETSc manual-page search; RAG (needs NVIDIA key + rag-data). |
| `agent_cards/` | | Only 4 cards: code_generator, compile_run, documentation, rag. |

### 0.2 `petsc_mcp_servers_tokamak` (github.com/engineer-scientist/..., cloned at /workspace/research/petsc_mcp_servers_tokamak, HEAD 8265e9c)
Achieved (Jul 23 – Aug 7 2026, ANL CELS GCE node, Claude Opus via ANL Argo):
* **Grad–Shafranov (steady, elliptic) only.** Agent-generated DMDA+SNES solver (267 lines), MMS-verified **p = 2.00** (max-norm 8.4e-4 @33² -> 1.3e-5 @257²), runs on 1 and 4 ranks (`artifacts/run-20260723-113024/verification.json`).
* Built-in orchestrator end-to-end -> DMPLEX+PetscFE GS solver; found/fixed the `cntlimit` false-failure bug (`artifacts/orchestrator-20260724{,-fixed}`, `patches/0004`).
* Shaped Solov'ev (Cerfon–Freidberg) equilibria ITER / NSTX / X-point: p = 2.00 on all three; finest max-norm 1.1e-7 (ITER); 1-vs-4-rank rel. diff 7e-4; q-profile vs FreeGS < 0.2 % (`artifacts/run-20260805-155724/`, metrics.md: 252 agent lines, 0 hand-written, 311 s wall, 25 LLM loops / 6 tool calls).
* 4 upstream patches, 9 agent design docs (incl. Visualization & Analysis agent — design only, no server), poster (US-RSE'26) + ANL intern slides (13 slides, 2026-07-29).
* ROADMAP milestone 10 (time-dependent / resistive MHD) **not started** in the repo.

### 0.3 Prior work found on ALCF (read-only inspection)
* **Polaris `~/NERCCS_2026/`** (Mar 2026): `plasma_symplectic.c` (86-line TSDISCGRAD-style harmonic-oscillator toy: H = ½|u|², block-symplectic S) and PBS weak-scaling scripts running **PETSc's own `src/ts/tutorials/hamiltonian/ex2.c`** (DMSwarm Vlasov–Poisson two-stream, `-ts_type basicsymplectic`) on 1–4 nodes, `preemptable` queue, PETSc 3.22.0 from `/soft/spack`. Two-stream energy/E-field logs exist; `logs/` (scaling outputs) is **empty** -> no scaling data recovered.
* **Aurora `~/mhd_project/mhd_scale.c`** (Apr 2026): DMPlex+TS skeleton whose IFunction body is empty ("GPU Kernel Launch Point") — not real physics.
* **Aurora `~/grokbot-tokamak/`, `~/grokbot-hello/`** (Oct 10–11 2026, from an earlier Grok Bot session): Python (non-PETSc) reduced-MHD study, 23 tasks on 1 node `debug` (job 8919984), report PDF + animation; SYCL hello-world. Useful physics context but not PETSc, not structure-preserving, not scaled.

### 0.4 Gap vs. the abstract
| Abstract claim | Status |
|---|---|
| Symplectic integrators in PETSc conserving phase-space volume / energy / momentum | Only a toy (NERCCS) + PETSc tutorial ex2 runs; **no own verified code, no long-time drift study vs RK** |
| Divergence-free field solver (∇·B = 0 discretely) | **Nothing exists** |
| Orchestrator + 4 specialists (modeling, NA, HPC codegen+exec, **visualization/analysis**) | 3 of 4 exist upstream; **viz/analysis agent design-only**; NA agent has no structure-preserving knowledge |
| Verification embedded in generated code (conservation checks, symplecticity tests) | Only MMS order-of-accuracy for GS. **No conservation/symplecticity checks** |
| Weak scaling on Polaris + Aurora, "thousands of GPU nodes" | **No scaling data at all.** Note Polaris has only ~560 nodes (max queue 496) — "thousands" is only possible on Aurora (large queue >= 2000 nodes) |

---------------------------------------------------------------------------------------------------
## 1. ALCF environment findings (read-only, 2026-10-11 ~20:40 ET)

### Aurora (`aurora-uan-0010`, user `$USER`)
* Default env: gcc/14.3.0, oneapi/release/2026.1.0, mpich 5.0.0 (aurora_test), cray-pals.
* **PETSc modules: `petsc/3.25.2-sycl` (default; built with Kokkos 5.1.1-sycl + kokkos-kernels, HDF5, CGNS, METIS/ParMETIS)** and `petsc/3.24.5-cpu`. Also hypre-sycl, ginkgo-sycl, amrex/warpx 26.08, cabana-sycl, kokkos (via petsc deps), heffte, superlu-dist.
* User build `~/petsc-aurora` (PETSc 3.25.0, `arch-linux-c-opt`, `--with-sycl=1 --with-syclc=icpx --with-debugging=0 --with-fc=0`, built 2026-04-30; libpetsc.so present). No Kokkos in that build -> prefer the module for GPU (Kokkos/SYCL) back ends.
* Project dir `/lus/flare/projects/DD_plasma_physics` exists but is **empty** and `myquota` shows **quota 1M, grace "expired"** -> effectively no project storage. Home: 2.4 G / 150 G used.
* Queues (PBS limits): `debug` 1–2 nodes/1 h; **`debug-scaling` 1–256 nodes/1 h**; `prod` routes to `small` 256–1024 nodes/12 h, `medium` 1025–1999/18 h, `large` 2000–10624/24 h; `capacity` 1–16 nodes/168 h. Note: prod minimum is 256 nodes, so 1–256 node points must go through debug/debug-scaling/capacity.
* Node: 6× Intel Data Center GPU Max (12 tiles), 2× SPR CPUs (public spec; to confirm with `-l select` test).

### Polaris (`polaris-login-02`)
* Default env: PrgEnv-nvidia 8.7.0 (nvidia/25.9), cray-mpich 9.1.0, cuda/12.9 & 13.0 available, craype-accel-nvidia80, Kokkos 4.3.01 / 4.6.02 in `/soft/libraries/kokkos`, spack-pe-* module trees.
* **No `petsc` module in the default view.** Usable PETSc installs:
  - `/soft/spack/hpc/0.8.1/.../petsc-3.22.0-5wsnxe4...` (used by NERCCS scripts; old: 3.22).
  - `~/petsc-hackathon/arch-linux-c-debug` (PETSc 3.25.0, `--with-cuda=1 --with-cuda-arch=80 --with-debugging=1`) — **debug build, unsuitable for timing.**
  - `~/petsc/arch-polaris-cuda-opt` — **configure FAILED 2026-08-07** ("CMake not found, needed to build METIS; rerun with --download-cmake"); no libpetsc.
  -> Need a fresh optimized CUDA(+Kokkos) PETSc 3.25.x/3.26 build on Polaris (fix: `module load cmake`/`--download-cmake`).
* Project dir `/eagle/DD_plasma_physics` (=/eagle/projects/...): contains `CSE674_MHD_project` (diffusion-model MHD ML work) and empty `$USER/tokamak-hpc`. **`myquota`: 75.83 G used vs 1M quota (over)** -> writes to /eagle likely blocked. Home: 11.6 G / 45 G.
* Queues: `debug` 1–2 nodes/1 h; `debug-scaling` 1–10/1 h; `prod` -> `small` 10–24/3 h, `medium` 25–99/6 h, `large` 100–496/24 h; `preemptable` 1–10/72 h; `capacity` 1–4/168 h.
* Node: 4× A100-40GB, 1× AMD Milan (public spec).

### Action items needing Sarthak (not done — external/irreversible)
1. **Storage**: ask ALCF support (support@alcf.anl.gov) to set project quotas on /lus/flare and /lus/eagle for DD_plasma_physics (both currently show 1M). Until then run from home dirs (fine for scaling runs that write only logs).
2. Approve any job submission / file writes on ALCF (Phase 2+).
3. Confirm which LLM backend the agents may use outside ANL CELS (Argo is ANL-internal; on this box we'd need an Anthropic key or equivalent).

---------------------------------------------------------------------------------------------------
## 2. (a) Structure-preserving methods to implement and analyze

Two application codes ("SP-PIC" and "SP-MHD") plus a particle-pusher micro-benchmark. All in C on PETSc, GPU-capable through PETSc Vec/Mat back ends (Kokkos on both machines; CUDA on Polaris, SYCL on Aurora).

### A1. Particle pushers (micro-benchmark, verification-heavy) — DMSwarm, no field solve
Single particle / ensemble in prescribed fields (uniform B; tokamak-like B from the Solov'ev ψ already in the tokamak repo; E×B; magnetic mirror).
* **Boris** (volume-preserving, not symplectic; gold standard)
* **Implicit midpoint / TSDISCGRAD** (energy-conserving discrete gradient; `u_t = S(u)∇F(u)` form)
* **TSBASICSYMPLECTIC** (Störmer–Verlet `-ts_basicsymplectic_type 2`, 3rd/4th order types) for separable Hamiltonians (electrostatic)
* Explicit/symplectic splitting for magnetized motion (e.g. He et al. explicit symplectic, Hamiltonian splitting) — optional
* Baseline: **TSRK 4** (RK4) and TSRK 5dp adaptive
Verification metrics: relative energy drift |H(t)−H(0)|/H(0) over 10^6–10^7 gyro-periods (bounded oscillation vs secular growth); canonical/toroidal angular momentum p_φ conservation in axisymmetric field (Noether); phase-space volume (Jacobian det of one-step map via finite differences / Liouville test on a particle cloud); gyro-phase error; convergence order. Expected: RK4 drift grows ~linearly in t; Boris/DG bounded.

### A2. SP-PIC: electrostatic Vlasov–Poisson (1D-1V -> 2D-2V) — DMSwarm + DMPlex FEM
Start from PETSc `src/ts/tutorials/hamiltonian/ex2.c` (Landau damping / two-stream, already run by Sarthak on Polaris) and ex3/ex4; build own version with:
* `TSBASICSYMPLECTIC` vs `TSDISCGRAD` vs RK4 for the particle ODE coupled to Poisson (`em_type primal`/mixed).
* Diagnostics: total energy (kinetic + field) drift, momentum, charge conservation, Landau damping rate vs linear theory (γ ≈ −0.1533 for k=0.5), two-stream growth rate.
* GPU: particle push via Kokkos (DMSwarm fields on device, PETSc 3.25 has device support partially — **verify** in 3.25.2 module), Poisson via PCGAMG on GPU.
* This is the **primary weak-scaling app** (particles/rank and cells/rank fixed).

### A3. SP-MHD / Maxwell: ∇·B = 0 to machine precision
* **Constrained transport (Yee/staggered) on DMSTAG**: B on faces, E on edges, B^{n+1} = B^n − Δt curl_h E; discrete div∘curl = 0 identity -> ∇_h·B stays at round-off (~1e-14) for all time. Test problems: 2D/3D Orszag–Tang (ideal MHD, with a simple FV/FD fluid part) or, simpler and more cleanly structure-preserving, **Maxwell / resistive induction equation with prescribed velocity** (field-loop advection, tokamak-geometry induction).
* **FEEC alternative on DMPlex**: H(curl)/H(div) Nédélec/Raviart–Thomas pair (PetscFE supports these on simplices — **verify element availability in 3.25**). Implicit-midpoint time stepping conserves discrete EM energy for Maxwell.
* Comparison baseline: cell-centered collocated FD for B (∇·B grows), with/without divergence cleaning.
* Verification: max|∇_h·B| vs time (target ≤ 1e-12 relative), EM energy conservation (implicit midpoint: exact up to solver tolerance; leapfrog: bounded), MMS convergence order for B and E.
* Secondary weak-scaling app (structured grid, cleanest scaling).

### A4. Analysis deliverables
Long-time drift plots (log-log energy error vs t, 5 methods), backward-error/modified-Hamiltonian argument for symplectic bounded drift, cost-vs-accuracy (energy error at fixed wall time — the fair comparison vs RK4), div B history, order-of-accuracy tables. Each is a slide.

---------------------------------------------------------------------------------------------------
## 3. (b) Tie-in to the multi-agent system

1. **NA agent upgrade** (`na_mcp_server.py`): new tool `select_structure_preserving_scheme(spec)` returning invariants to preserve (energy, momentum, phase-space volume, div B) and a deterministic map invariant -> PETSc choice (TSBASICSYMPLECTIC / TSDISCGRAD / Boris kernel; DMSTAG CT / FEEC Nédélec–RT; DMSWARM). Today the NA prompt only chooses grid/discretization/integrator.
2. **Verification-contract in code-gen**: the code-gen query must require the generated program to print a machine-readable `VERIFY {...}` JSON line (energy drift, |div B|max, momentum error, MMS errors). Compile-run checks it, not just compilation.
3. **New Visualization & Analysis agent** (design exists in `design/Visualization_and_Analysis_Agent_Design_Document.md`, port 8087): parses VERIFY lines + PETSc `-log_view`; produces drift / divB / scaling plots; emits pass/fail + rollback signal to the orchestrator. Agent card `agent_cards/visualization_analysis.md`.
4. **HPC execution extension** (compile-run): optional PBS job-script generator for Polaris/Aurora (template from this plan's scripts); submission stays **human-approved**.
5. New agent cards: `numerical_analysis.md` (with SP skills), `visualization_analysis.md`, `orchestrator.md`, `pde_modeling.md` (repo currently has only 4 cards).
6. **Experiment**: run the pipeline on 4 prompts (Boris pusher in tokamak field; Vlasov–Poisson Landau damping with symplectic TS; CT induction equation; Maxwell FEEC). Metrics as in the tokamak repo: pass rate, LLM loops, tool calls, human-edited lines, verification outcome (incl. "agent produced RK4 -> analysis agent flags drift -> rollback -> symplectic" demo). Hand-written reference codes (A1–A3) are the oracle.
Constraint: LLM backend (Argo is ANL-internal) — flagged unknown.

---------------------------------------------------------------------------------------------------
## 4. (c) Weak-scaling study design

Apps: **SP-PIC** (primary; fixed particles/GPU, e.g. 2–8 M particles/GPU-tile and fixed cells/GPU) and **SP-MHD CT** (fixed DMSTAG cells/GPU, e.g. 256³/GPU on A100, tune for PVC tile). Measure time per step (median of N steps after warm-up), per-stage `-log_view` (push, deposit, Poisson KSP/PCGAMG, halo/migration), GAMG iterations, efficiency E(N) = T(1)/T(N). Ranks: 1 per GPU (Polaris 4/node), 1 per tile (Aurora 12/node) — confirm with a 1-node test. Also a "CPU-only" point per machine for context.

Job length: ~15 min wall per point (incl. startup), single run per point after a 1–2 node pilot; 3 repetitions only for ≤ 64 nodes (cheap) to show variance.

### Polaris (4,964 node-h left)
Nodes 1, 2, 4, 8 (debug/debug-scaling) · 16, 32, 64 (small/medium) · 128, 256 (large) [optional 448].
Per app per sweep: Σ nodes = 511 × 0.25 h ≈ 128 node-h. 2 apps × (1 sweep + small-scale repeats) ≈ 330 node-h. Dev/pilot/profiling ≈ 250 node-h. Optional 448-node point ×2 apps ≈ 225. **Polaris total ≈ 600–800 node-h (≤ 16 % of remaining).**

### Aurora (11,999 node-h left)
Nodes 1, 2 (debug) · 4…256 (debug-scaling, 1 h cap) · 512, 1024 (small) · optional 2048 (large; this is the only "thousands of GPU nodes" point).
Per app per sweep to 1024: Σ = 2047 × 0.25 h ≈ 512 node-h. 2 apps ≈ 1,025; repeats ≤ 64 nodes ≈ 50; dev/pilot (SYCL/Kokkos bring-up, 1–8 nodes) ≈ 400; one 2048-node point per app ≈ 2 × 2048 × 0.25 ≈ 1,025. **Aurora total ≈ 1,500 (to 1024 nodes) – 2,600 (incl. 2048) node-h (≤ 22 % of remaining).** Buffer kept for failures/re-runs. Risk: large-queue wait times before Nov 7 — submit 1024/2048 by ~Oct 30.

Weak-scaling deliverable: efficiency plot (both machines, both apps), stage breakdown bars, conservation metrics shown to be **independent of node count** (energy drift / div B unchanged at scale — key structure-preservation-at-scale message).

---------------------------------------------------------------------------------------------------
## 5. (d) Milestones (dates ET)

| # | Dates | Milestone | Where |
|---|---|---|---|
| M0 | Oct 11–13 | This survey/plan; Sarthak approves plan, files storage-quota request; decide LLM backend | — |
| M1 | Oct 12–16 | Build PETSc (CPU) on Grok Bot box; A1 pusher suite + drift/symplecticity tests; reproduce ex2 Landau damping rate | box |
| M2 | Oct 14–20 | A3 DMSTAG CT induction/Maxwell: div B ~1e-14, energy, MMS; A2 own SP-PIC with TS comparison | box |
| M3 | Oct 17–22 | Polaris: build optimized CUDA(+Kokkos) PETSc 3.25.x; Aurora: validate `petsc/3.25.2-sycl` module; 1–2 node GPU pilots (needs approval: writes + debug jobs) | ALCF |
| M4 | Oct 20–27 | Agent work: NA SP tool, VERIFY contract, viz/analysis server, 4 agent cards; run 4-prompt experiment, collect metrics | box |
| M5 | Oct 23–Nov 2 | Weak-scaling sweeps: Polaris 1–256, Aurora 1–1024 (+2048 submit by Oct 30) | ALCF |
| M6 | Nov 2–7 | Buffer / re-runs; archive all data to home + GitHub; allocation ends Nov 7 | ALCF |
| M7 | Nov 8 → talk | Analysis write-up; 12-min CCS 2026 deck (~10–12 slides): motivation (complex system) · why SP methods (drift plot) · div B · agent hierarchy + verification loop · agent experiment metrics · weak scaling both machines · conclusions | box |

Talk outline (12 min ≈ 11 slides): 1 title · 2 fusion plasma as complex system · 3 failure of RK over long times · 4 symplectic/DG pushers result · 5 CT/FEEC div B result · 6 SP-PIC Landau damping · 7 multi-agent architecture (from tokamak repo figure) · 8 verification-driven generation + metrics (incl. GS p=2.00 prior result) · 9 weak scaling Polaris · 10 weak scaling Aurora · 11 conclusions/next.

---------------------------------------------------------------------------------------------------
## 6. Unknowns / risks
* Exact GPU support of DMSwarm particle operations in PETSc 3.25.2 (Kokkos/SYCL) — may need custom Kokkos kernels for the push; verify early (M3).
* PetscFE Nédélec/RT availability and robustness for 3D FEEC in 3.25 — CT on DMSTAG is the safe default.
* Project storage quotas (both 1M) — must be fixed or stay in home dirs.
* Aurora large-queue (≥2000 nodes) access/wait time for a DD project before Nov 7; "thousands of GPU nodes" claim may need softening to "up to 1024 Aurora nodes (12,288 GPU tiles)".
* LLM backend for agents off ANL network.
* Box resources: 8 cores, 15 GB RAM (shared with other agents) — fine for development, not for large runs.
* CCS 2026 talk date not yet known.
