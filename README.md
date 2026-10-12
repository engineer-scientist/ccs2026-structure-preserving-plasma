# Structure-preserving plasma simulation + PETSc multi-agent AI (CCS 2026)

Companion repository for the talk **"Generative AI for Scalable and Structure-Preserving Plasma Simulations: A Complex Systems Approach to Nuclear Fusion Energy"**, Sarthak Sharma and Matthew Knepley, University at Buffalo. The talk is virtual, at the Conference on Complex Systems 2026, on 13 Oct 2026.

It holds the code, data, figures, slides and report for a short study of *structure-preserving* numerical methods in PETSc (methods that keep the physical invariants a plasma code must not lose) and of how they connect to an AI agent pipeline that generates PETSc codes.

## Layout
| path | contents |
|---|---|
| `presentation/` | 14-slide, 12-minute talk deck (`.pptx`, incl. strong scaling and long-run conservation; Takeaways last), with speaker notes |
| `report/` | Full project report (PDF, 46 pp., incl. Phase 4 chapter) written for non-experts |
| `docs/PLAN.md` | Original research plan and survey |
| `results/RESULTS.md` | **Main results write-up** (phases 2 and 3), with numbers and caveats |
| `results/figs/` | All figures used in the talk and report |
| `results/pusher/`, `results/*.csv`, `results/*.json` | Raw particle-pusher, constrained-transport, Landau and scaling data |
| `results/alcf/` | ALCF pilot job summaries, plus PBS and build scripts (`results/alcf/scripts/`) |
| `code/sp/pusher/` | Single-particle pushers in PETSc TS (RK4, TSBASICSYMPLECTIC, Boris, TSDISCGRAD) |
| `code/sp/ct/`, `code/sp/ct_stag/` | Constrained-transport induction: numpy and PETSc DMSTAG versions |
| `code/sp/ex2/` | PETSc `ex2` (DMSwarm Vlasov–Poisson, Landau damping) run configuration and analysis |
| `code/verify_contract/` | `VERIFY{}` conservation-contract prototype and patches for the agent pipeline |
| `phase3/` | Load-imbalance diagnosis, `ex2` grid-hash patch (`patch/ex2_gb_wrap_hash.patch`), PBS scripts, weak-scaling data and logs, `NOTES.md` |

Paths inside `RESULTS.md` such as `sp/pusher/pusher.c` refer to `code/sp/...` in this repository.

## Key results (details in `results/RESULTS.md`)
- **Particle pushers (10^6 steps).** RK4 drifts secularly: the energy error reaches 19% in a Solov'ev tokamak field and 52% in uniform B. Boris and TSDISCGRAD keep the invariants bounded or at round-off (for example, Boris energy error 1.3e-13 in the Solov'ev field).
- **TSBASICSYMPLECTIC caveat.** PETSc's TSBASICSYMPLECTIC assumes a separable Hamiltonian. With a v×B force the Hamiltonian is not separable, and the method blows up in fewer than 2,000 steps. It is correct for electrostatic Vlasov–Poisson, which is the case `ex2` uses.
- **∇·B = 0.** Constrained transport on DMSTAG keeps div B at round-off (≤ 3e-24). A cell-centred B advected with the same scheme reaches about 5e-2.
- **Landau damping on ALCF Aurora and Polaris.** The fitted damping rate is within about 5% of linear theory (−0.1533).
- **Scaling fix.** In `DMPlexLocatePoints`, the rank that owns the periodic wrap cell was 18–41× slower than the others. A grid-hash patch (`-gb_wrap_hash_factor 2`) fixes this and makes runs **12.2× faster per step at 16 Aurora nodes**. Weak-scaling efficiency is 75% from 1 to 16 nodes on Aurora and 80% from 1 to 8 nodes on Polaris. Particle count and charge are conserved exactly at every scale.
- **VERIFY{} contract.** Generated codes must print a JSON line of invariants, and the pipeline loops until the tolerances pass. TSBASICSYMPLECTIC in a B field correctly FAILs.

## Reproducing
1. Build PETSc 3.26.0 (`./configure --with-debugging=0`, plus `--with-cuda` / SYCL on Polaris / Aurora; see `results/alcf/scripts/build_petsc_aurora.sh` and `phase3/scripts/rebuild_polaris.sh`).
2. Pushers: `cd code/sp/pusher && make PETSC_DIR=... && ./run_all.sh && python plot_pusher.py`.
3. Constrained transport: `python code/sp/ct/ct2d.py`; `cd code/sp/ct_stag && make && ./ct_stag`.
4. Landau / ex2: use PETSc `src/ts/tutorials/hamiltonian/ex2.c` with the options in `code/sp/ex2/*.opts` and `common_opts.txt`. For the scaling fix, apply `phase3/patch/ex2_gb_wrap_hash.patch`.
5. ALCF: the PBS scripts are in `results/alcf/scripts/` and `phase3/scripts/` (`ws_aurora*.pbs`, `ws_polaris*.pbs`). Adjust the project, queue and paths to your environment. Reduce the logs with `phase3/scripts/reduce3.py` and plot them with `make_fig.py`.

User names and home paths in the logs and scripts have been replaced with `$USER` / `$HOME`.

## Related repositories
- PETSc MCP multi-agent servers: https://gitlab.com/petsc/petsc_mcp_servers
- Tokamak application: https://github.com/engineer-scientist/petsc_mcp_servers_tokamak

## Acknowledgement
This research used resources of the Argonne Leadership Computing Facility (ALCF), a U.S. DOE Office of Science user facility at Argonne National Laboratory, through the allocation **DD_plasma_physics** (Aurora and Polaris).

## License
No license has been chosen yet, so all rights are reserved by default. Files derived from PETSc (`code/sp/ex2/ex2.c`, `phase3/patch/ex2.c`) remain under PETSc's BSD-2-Clause license.

## Phase 4 (Oct 12, 2026): strong scaling, 64-node weak scaling, long-run conservation
See `phase4/` and the Phase 4 section of `results/RESULTS.md`.
- Strong scaling (9.83M particles): Aurora 1→16 nodes 88.6→19.1 s (4.6×, 29 % eff.); Polaris 1→8 nodes 235→43.7 s (5.4×, 67 %).
- Weak scaling vs 1 node (TSStep): 63 % at 32 nodes, 58 % at 64 nodes (6,144 ranks, 315M particles) on Aurora.
- Long Landau runs: energy error bounded (~4e-6) at Δt=0.01 to t=60, but secular numerical heating at Δt=0.05 in long runs (0.44 % by t=300; 81 % for a 3,200-particle run to t=3000).
- Pushers, 10⁸ steps: Boris bounded (Solov'ev ΔH ~1e-12), RK4 drifts 28–67 %; TSDISCGRAD stalled (not a 10⁸-step result).
